package main

import "os"

type BitWriter struct {
	file                   *os.File
	word, word_bits        uint32
	buffer                 []byte
	buffer_pos, buffer_len int
}

func BitWriterOpen(filename string) (*BitWriter, error) {
	if f, err := os.Create(filename); err != nil {
		return nil, err
	} else {
		return &BitWriter{file: f, buffer: make([]byte, 1<<12), buffer_len: 1 << 12}, nil
	}
}

func (br *BitWriter) Close(flush bool) error {
	if br.file != nil {
		defer (func() {
			br.file = nil
		})()

		if flush {
			for i := 0; i < 4; i++ {
				br.WriteBits(0, 8)
			}
		}

		if br.buffer_pos > 0 {
			br.file.Write(br.buffer[:br.buffer_pos])
		}

		return br.file.Close()
	} else {
		return nil
	}
}

func (br *BitWriter) refill() {
	for br.word_bits >= 8 {
		if br.buffer_pos == br.buffer_len {
			_, err := br.file.Write(br.buffer[:br.buffer_pos])
			br.buffer_pos = 0

			if err != nil {
				panic(err)
			}
		}

		br.buffer[br.buffer_pos] = uint8(br.word >> 24)
		br.word <<= 8
		br.buffer_pos++
		br.word_bits -= 8
	}
}

func (br *BitWriter) WriteBits(v, nb uint32) {
	br.word |= v << (32 - br.word_bits - nb)
	br.word_bits += nb

	br.refill()
}
