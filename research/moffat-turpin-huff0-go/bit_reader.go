package main

import "os"

type BitReader struct {
	file                   *os.File
	word, word_bits        uint32
	buffer                 []byte
	buffer_pos, buffer_len int
}

func BitReaderOpen(filename string) (*BitReader, error) {
	if f, err := os.Open(filename); err != nil {
		return nil, err
	} else {
		return &BitReader{file: f, buffer: make([]byte, 1<<12)}, nil
	}
}

func (br *BitReader) Close() error {
	if br.file != nil {
		defer (func() {
			br.file = nil
		})()

		return br.file.Close()
	} else {
		return nil
	}
}

func (br *BitReader) refill() {
	for br.word_bits < 24 {
		if br.buffer_pos == br.buffer_len {
			n, err := br.file.Read(br.buffer)

			br.buffer_pos = 0
			br.buffer_len = n
			if n == 0 {
				return
			}

			if err != nil {
				panic(err)
			}
		}

		br.word |= uint32(br.buffer[br.buffer_pos]) << (24 - br.word_bits)
		br.buffer_pos++
		br.word_bits += 8
	}
}

func (br *BitReader) PeekBits(nb uint32) uint32 {
	br.refill()
	return br.word >> (32 - nb)
}

func (br *BitReader) Discard(nb uint32) {
	if nb > br.word_bits {
		br.word = 0
		br.word_bits = 0
	} else {
		br.word <<= nb
		br.word_bits -= nb
	}
}

func (br *BitReader) EOF() bool {
	br.refill()
	return br.word_bits == 0
}
