use std::{
  fs::File,
  io::{self, Read, Write},
};

const BIT_BUF_SIZE: u32 = 1 << 12;

pub struct BitReader<'f> {
  file: &'f File,
  buffer: [u8; BIT_BUF_SIZE as usize],
  buffer_pos: u32,
  buffer_len: u32,
  word: u32,
  word_bits: u32,
}

impl<'f> BitReader<'f> {
  pub fn new(file: &'f File) -> Self {
    BitReader {
      file,
      buffer: [0; BIT_BUF_SIZE as usize],
      buffer_pos: 0,
      buffer_len: 0,
      word: 0,
      word_bits: 0,
    }
  }

  fn refill(&mut self) -> io::Result<usize> {
    let mut bytes_read: usize = 0;

    while self.word_bits < 24 {
      if self.buffer_pos == self.buffer_len {
        self.buffer_len = self.file.read(&mut self.buffer)? as u32;
        self.buffer_pos = 0;
        if self.buffer_pos == self.buffer_len {
          break;
        }
      }

      bytes_read += 1;

      self.word |= (self.buffer[self.buffer_pos as usize] as u32) << (24 - self.word_bits);
      self.buffer_pos += 1;
      self.word_bits += 8;
    }

    Ok(bytes_read)
  }

  pub fn peek(&mut self, nb: u32) -> io::Result<u32> {
    self.refill()?;
    Ok(self.word >> (32 - nb))
  }

  pub fn discard(&mut self, nb: u32) -> () {
    if nb < self.word_bits {
      self.word_bits -= nb;
      self.word <<= nb;
    } else {
      self.word_bits = 0;
      self.word = 0;
    }
  }

  pub fn read_bits(&mut self, nb: u32) -> io::Result<u32> {
    let v = self.peek(nb)?;
    self.discard(nb);

    Ok(v)
  }
}

pub struct BitWriter<'f> {
  file: &'f File,
  buffer: [u8; BIT_BUF_SIZE as usize],
  buffer_pos: u32,
  word: u32,
  word_bits: u32,
}

impl<'f> BitWriter<'f> {
  pub fn new(file: &'f File) -> Self {
    BitWriter {
      file,
      buffer: [0; BIT_BUF_SIZE as usize],
      buffer_pos: 0,
      word: 0,
      word_bits: 0,
    }
  }

  fn refill(&mut self) -> io::Result<usize> {
    let mut bytes_written: usize = 0;

    while self.word_bits >= 8 {
      if self.buffer_pos == BIT_BUF_SIZE {
        self.file.write(&self.buffer[0..self.buffer_pos as usize])?;
        self.buffer_pos = 0;
      }

      bytes_written += 1;

      self.buffer[self.buffer_pos as usize] = (self.word >> 24) as u8;
      self.word <<= 8;
      self.buffer_pos += 1;
      self.word_bits -= 8;
    }

    Ok(bytes_written)
  }

  pub fn write_bits(&mut self, v: u32, nb: u32) -> io::Result<usize> {
    assert!(v < (1 << nb));
    self.word |= v << (32 - self.word_bits - nb);
    self.word_bits += nb;
    self.refill()?;

    Ok(nb as usize)
  }

  pub fn flush(&mut self) -> io::Result<()> {
    for _ in 0..4 {
      self.write_bits(0, 8)?;
    }

    if self.buffer_pos > 0 {
      self.file.write(&self.buffer[0..self.buffer_pos as usize])?;
      self.buffer_pos = 0;
    }

    Ok(())
  }
}
