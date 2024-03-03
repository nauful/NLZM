package main

import (
	"flag"
	"fmt"
	"os"
)

const (
	INITIAL_FRAME_SIZE = 1 << 10
	FRAME_MAX_SIZE     = 1 << 15
)

func usage() {
	fmt.Fprintf(os.Stderr, "usage: main c/d [input] [output]\n")
	flag.PrintDefaults()
	os.Exit(2)
}

func main() {
	// compress_file("test.txt", "out.bin")
	// decompress_file("out.bin", "out.txt")

	flag.Usage = usage
	flag.Parse()

	args := flag.Args()
	if len(args) < 3 || (args[0] != "c" && args[0] != "d") {
		usage()
		os.Exit(1)
	}

	switch args[0] {
	case "c":
		compress_file(args[1], args[2])
	case "d":
		decompress_file(args[1], args[2])
	}
}

func compress_file(input_filename, output_filename string) {
	fin, fin_err := BitReaderOpen(input_filename)
	if fin_err != nil {
		panic(fin_err)
	}
	defer fin.Close()

	fout, fout_err := BitWriterOpen(output_filename)
	if fout_err != nil {
		panic(fout_err)
	}
	defer fout.Close(true)

	frame := make([]uint8, FRAME_MAX_SIZE)
	frame_limit := uint32(INITIAL_FRAME_SIZE)
	frame_pos := uint32(0)
	code_table := CreateCodeTable(frame[0:0])

	tot_nb := uint32(0)
	for !fin.EOF() {
		frame[frame_pos] = uint8(fin.PeekBits(8))
		frame_pos++
		fin.Discard(8)

		if frame_pos == frame_limit {
			fout.WriteBits(frame_pos&0xFF, 8)
			fout.WriteBits((frame_pos>>8)&0xFF, 8)
			for i := uint32(0); i < frame_pos; i++ {
				nb := code_table.code_len[frame[i]]
				tot_nb += nb
				fout.WriteBits(code_table.code[frame[i]], nb)
			}

			code_table = CreateCodeTable(frame[:frame_pos])

			frame_pos = 0
			if frame_limit < FRAME_MAX_SIZE {
				frame_limit <<= 1
			}
		}
	}

	fout.WriteBits(frame_pos&0xFF, 8)
	fout.WriteBits((frame_pos>>8)&0xFF, 8)
	for i := uint32(0); i < frame_pos; i++ {
		nb := code_table.code_len[frame[i]]
		tot_nb += nb
		fout.WriteBits(code_table.code[frame[i]], nb)
	}

	if frame_pos > 0 {
		frame_pos = 0
		fout.WriteBits(frame_pos&0xFF, 8)
		fout.WriteBits((frame_pos>>8)&0xFF, 8)
	}
}

func decompress_file(input_filename, output_filename string) {
	fin, fin_err := BitReaderOpen(input_filename)
	if fin_err != nil {
		panic(fin_err)
	}
	defer fin.Close()

	fout, fout_err := BitWriterOpen(output_filename)
	if fout_err != nil {
		panic(fout_err)
	}

	defer fout.Close(false)

	frame := make([]uint8, FRAME_MAX_SIZE)
	frame_pos := uint32(0)
	decode_table := CreateDecodeTable(frame[0:0])

	for !fin.EOF() {
		frame_pos = fin.PeekBits(8)
		fin.Discard(8)
		frame_pos |= fin.PeekBits(8) << 8
		fin.Discard(8)

		if frame_pos == 0 {
			break
		}

		for i := uint32(0); i < frame_pos; i++ {
			word := fin.PeekBits(CODE_LENGTH_LIMIT)
			word_len := uint32(1)
			for word >= decode_table.base_code[word_len+1] {
				word_len++
			}
			word_idx := (word - decode_table.base_code[word_len]) >> (CODE_LENGTH_LIMIT - word_len)
			sym := decode_table.sym_table[word_idx+decode_table.base_offset[word_len]]

			frame[i] = uint8(sym)
			fout.WriteBits(uint32(sym), 8)
			fin.Discard(word_len)
		}

		decode_table = CreateDecodeTable(frame[:frame_pos])
	}
}
