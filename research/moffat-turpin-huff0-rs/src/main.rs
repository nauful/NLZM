use std::{
  env,
  fs::{File, OpenOptions},
  io::{self, Read, Seek, Write},
  path::Path,
  time::Instant,
};

mod bit_io;
mod coder;

fn main() {
  let args: Vec<String> = env::args().collect();
  if args.len() != 4 || args[1] != "c" || args[1] != "d" {
    eprintln!("Usage: {} c/d input output", args[0]);
    std::process::exit(1);
  }

  match args[1].as_str() {
      "c" => compress_file(Path::new(&args[2]), Path::new(&args[3]))
        .expect("Failed to compress the file"),
      "d" => decompress_file(Path::new(&args[2]), Path::new(&args[3]))
        .expect("Failed to dempress the file"),
      _ => unreachable!(),
  }
}

fn compress_file(input_path: &Path, out_path: &Path) -> io::Result<()> {
  let mut fin = File::open(input_path)?;
  let mut fout = OpenOptions::new()
    .write(true)
    .create(true)
    .truncate(true)
    .open(out_path)?;

  let mut frame = [0 as u8; coder::FRAME_MAX_SIZE as usize];
  let mut frame_limit: u32 = coder::FRAME_INITIAL_SIZE;
  let mut tot_nb: u32 = 0;
  let mut code_table: Vec<coder::CodingTableEntry> = coder::create_code_table(&frame[0..0]);
  let mut writer = bit_io::BitWriter::new(&fout);

  loop {
    let nb_read: usize = fin.read(&mut frame[0..frame_limit as usize])?;

    writer.write_bits(((nb_read >> 8) & 0xFF) as u32, 8)?;
    writer.write_bits((nb_read & 0xFF) as u32, 8)?;

    for i in 0..nb_read {
      tot_nb += code_table[frame[i as usize] as usize].code_len;
      writer.write_bits(
        code_table[frame[i] as usize].code,
        code_table[frame[i] as usize].code_len,
      )?;
    }

    if nb_read == 0 {
      break;
    }

    code_table = coder::create_code_table(&frame[0..nb_read]);
    if frame_limit < coder::FRAME_MAX_SIZE {
      frame_limit <<= 1;
    }
  }

  writer.flush()?;

  println!(
    "{}KB -> {}KB {} KB",
    fin.stream_position()? >> 10,
    fout.stream_position()? >> 10,
    tot_nb >> 13
  );

  return Ok(());
}

fn decompress_file(input_path: &Path, out_path: &Path) -> io::Result<()> {
  let mut fin = File::open(input_path)?;
  let mut fout = OpenOptions::new()
    .write(true)
    .create(true)
    .truncate(true)
    .open(out_path)?;

  let time0 = Instant::now();

  let mut frame = [0 as u8; coder::FRAME_MAX_SIZE as usize];
  let mut reader = bit_io::BitReader::new(&fin);
  let mut decode_table = coder::create_decode_table(&frame[0..0]);

  loop {
    let frame_limit: u32 = (reader.read_bits(8)? << 8) + reader.read_bits(8)?;
    if frame_limit == 0 {
      break;
    }

    for i in 0..frame_limit as usize {
      let word = reader.peek(coder::CODE_LENGTH_LIMIT)?;

      let mut word_len: u32 = 1;
      while word >= decode_table.base_code[(word_len + 1) as usize] {
        word_len += 1
      }

      let word_idx =
        (word - decode_table.base_code[word_len as usize]) >> (coder::CODE_LENGTH_LIMIT - word_len);
      let sym =
        decode_table.sym_table[(word_idx + decode_table.base_offset[word_len as usize]) as usize];

      frame[i] = sym as u8;
      reader.discard(word_len);
    }

    decode_table = coder::create_decode_table(&frame[0..frame_limit as usize]);

    fout.write(&frame[0..frame_limit as usize])?;
  }

  let time_elapsed = time0.elapsed();

  println!(
    "{}KB -> {}KB in {:?}",
    fin.stream_position()? >> 10,
    fout.stream_position()? >> 10,
    time_elapsed
  );

  Ok(())
}
