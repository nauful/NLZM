package main

import (
	"sort"
)

const (
	CODE_LENGTH_LIMIT = 14
	CODE_RANGE        = 0x100
)

type CodeEntry struct {
	sym, count, code_len uint32
}

type CodingTable struct {
	code, code_len []uint32
}

type DecodingTable struct {
	base_code, base_offset []uint32
	sym_table              []uint32
}

func CreateFreqTable(frame []uint8) []CodeEntry {
	entries := make([]CodeEntry, CODE_RANGE)

	for i := 0; i < CODE_RANGE; i++ {
		entries[i].sym = uint32(i)
		entries[i].count = 1
	}

	for i := 0; i < len(frame); i++ {
		entries[frame[i]].count++
	}

	count_compare := func(i, j int) bool { return entries[i].count < entries[j].count }
	sort.Slice(entries, count_compare)

	return entries
}

func AssignCodes(table CodingTable) {
	len_count := make([]uint32, 17)
	next_code := make([]uint32, 17)
	cur_code := uint32(0)

	max_len := uint32(0)
	for i := 0; i < CODE_RANGE; i++ {
		len_count[table.code_len[i]]++
		if table.code_len[i] > max_len {
			max_len = table.code_len[i]
		}
	}

	for bits := uint32(1); bits <= max_len; bits++ {
		cur_code += len_count[bits-1]
		cur_code <<= 1
		next_code[bits] = cur_code
	}
	cur_code += len_count[max_len]

	for i := 0; i < CODE_RANGE; i++ {
		table.code[i] = next_code[table.code_len[i]]
		next_code[table.code_len[i]]++
	}
}

func LeftAssignDecodes(table CodingTable, align_bits uint32) (base_code []uint32, base_offset []uint32) {
	len_count := make([]uint32, 17)
	base_code = make([]uint32, 17)
	base_offset = make([]uint32, 17)

	cur_code := uint32(0)

	max_len := uint32(0)
	for i := 0; i < CODE_RANGE; i++ {
		len_count[table.code_len[i]]++
		if table.code_len[i] > max_len {
			max_len = table.code_len[i]
		}
	}

	for bits := uint32(1); bits <= max_len; bits++ {
		base_offset[bits] = base_offset[bits-1] + len_count[bits-1]
		cur_code += len_count[bits-1]
		cur_code <<= 1

		base_code[bits] = cur_code << (align_bits - bits)
	}

	cur_code += len_count[max_len]
	base_code[max_len+1] = cur_code << (align_bits - max_len)

	for max_len < align_bits {
		base_code[max_len+2] = base_code[max_len+1]
		max_len++
	}

	return
}

func CreateCodeTable(frame []uint8) CodingTable {
	entries := CreateFreqTable(frame)

	tree_count := make([]uint32, 2*CODE_RANGE)
	bit_len := make([]uint32, 2*CODE_RANGE)
	left := make([]uint32, 2*CODE_RANGE)
	right := make([]uint32, 2*CODE_RANGE)

	for {
		p0 := uint32(0)
		p1 := uint32(CODE_RANGE + 1)

		total := uint32(0)
		for i := 0; i < CODE_RANGE; i++ {
			tree_count[i] = entries[i].count
			total += tree_count[i]
		}

		for i := CODE_RANGE; i < 2*CODE_RANGE; i++ {
			tree_count[i] = 0xFFFFFFFF
		}

		for w := CODE_RANGE + 1; w < 2*CODE_RANGE; w++ {
			if tree_count[p0] <= tree_count[p1] {
				left[w] = p0
				p0++
			} else {
				left[w] = p1
				p1++
			}

			if tree_count[p0] <= tree_count[p1] {
				right[w] = p0
				p0++
			} else {
				right[w] = p1
				p1++
			}

			tree_count[w] = tree_count[left[w]] + tree_count[right[w]]
		}

		bit_len[2*CODE_RANGE-1] = 0
		for i := 2*CODE_RANGE - 1; i > CODE_RANGE; i-- {
			bit_len[left[i]] = bit_len[i] + 1
			bit_len[right[i]] = bit_len[i] + 1
		}

		max_len := uint32(0)
		for i := 0; i < CODE_RANGE; i++ {
			if bit_len[i] > max_len {
				max_len = bit_len[i]
			}

			entries[i].code_len = bit_len[i]
		}

		if max_len <= CODE_LENGTH_LIMIT {
			break
		}

		bs := max_len - CODE_LENGTH_LIMIT
		for i := 0; i < CODE_RANGE; i++ {
			entries[i].count >>= bs
			if entries[i].count == 0 {
				entries[i].count = 1
			}
		}
	}

	table := CodingTable{
		code:     make([]uint32, CODE_RANGE),
		code_len: make([]uint32, CODE_RANGE),
	}

	for i := 0; i < CODE_RANGE; i++ {
		table.code_len[entries[i].sym] = entries[i].code_len
	}

	AssignCodes(table)
	return table
}

func CreateDecodeTable(frame []uint8) DecodingTable {
	table := CreateCodeTable(frame)

	base_code, base_offset := LeftAssignDecodes(table, CODE_LENGTH_LIMIT)

	sym_table := make([]uint32, CODE_RANGE)
	cur_offset := make([]uint32, cap(base_offset))
	copy(cur_offset, base_offset)
	for i := 0; i < CODE_RANGE; i++ {
		len := table.code_len[i]
		sym_table[cur_offset[len]] = uint32(i)
		cur_offset[len]++
	}

	return DecodingTable{
		base_code:   base_code,
		base_offset: base_offset,
		sym_table:   sym_table,
	}
}
