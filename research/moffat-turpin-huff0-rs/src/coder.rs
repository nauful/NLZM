pub const FRAME_INITIAL_SIZE: u32 = 1 << 12;
pub const FRAME_MAX_SIZE: u32 = 1 << 15;

pub const CODE_LENGTH_LIMIT: u32 = 14;
const CODE_RANGE: usize = 0x100;

#[derive(Debug)]
struct CodeEntry {
  sym: u32,
  count: u32,
  code_len: u32,
}

pub struct CodingTableEntry {
  pub code: u32,
  pub code_len: u32,
}

impl PartialEq for CodeEntry {
  fn eq(&self, other: &Self) -> bool {
    self.sym == other.sym
  }
}

impl Eq for CodeEntry {}

impl PartialOrd for CodeEntry {
  fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
    self.count.partial_cmp(&other.count)
  }
}

impl Ord for CodeEntry {
  fn cmp(&self, other: &Self) -> std::cmp::Ordering {
    self.count.cmp(&other.count)
  }
}

fn create_freq_table(frame: &[u8]) -> Vec<CodeEntry> {
  let mut entries: Vec<CodeEntry> = Vec::default();
  for i in 0..CODE_RANGE as u32 {
    entries.push(CodeEntry {
      sym: i,
      count: 1,
      code_len: 0,
    })
  }

  for c in frame {
    entries[*c as usize].count += 1;
  }

  entries.sort_by(|a, b| a.count.cmp(&b.count));

  entries
}

fn assign_codes(table: &mut Vec<CodingTableEntry>) {
  let mut len_count = [0; 17];
  let mut next_code = [0; 17];

  let mut max_len = 0;
  for code_len in table.iter().map(|e| e.code_len) {
    len_count[code_len as usize] += 1;
    max_len = u32::max(max_len, code_len);
  }

  let mut cur_code = 0;
  for bits in 1..=max_len {
    cur_code += len_count[(bits - 1) as usize];
    cur_code <<= 1;
    next_code[bits as usize] = cur_code;
  }

  for e in table.iter_mut() {
    e.code = next_code[e.code_len as usize];
    next_code[e.code_len as usize] += 1;
  }
}

fn left_assign_codes(table: &Vec<CodingTableEntry>, align_bits: u32) -> (Vec<u32>, Vec<u32>) {
  let mut len_count = [0; 17];
  let mut base_code = [0; 17];
  let mut base_offset = [0; 17];

  let mut max_len = 0;
  for code_len in table.iter().map(|e| e.code_len) {
    len_count[code_len as usize] += 1;
    max_len = u32::max(max_len, code_len);
  }

  let mut cur_code = 0;
  for bits in 1..=max_len {
    base_offset[bits as usize] = base_offset[(bits - 1) as usize] + len_count[(bits - 1) as usize];
    cur_code += len_count[(bits - 1) as usize];
    cur_code <<= 1;

    base_code[bits as usize] = cur_code << (align_bits - bits);
  }

  cur_code += len_count[max_len as usize];
  base_code[(max_len + 1) as usize] = cur_code << (align_bits - max_len);

  while max_len < align_bits {
    base_code[(max_len + 2) as usize] = base_code[(max_len + 1) as usize];
    max_len += 1;
  }

  (base_code.to_vec(), base_offset.to_vec())
}

pub fn create_code_table(frame: &[u8]) -> Vec<CodingTableEntry> {
  let mut entries = create_freq_table(frame);

  let mut tree_count = [0; 2 * CODE_RANGE];
  let mut bit_len = [0; 2 * CODE_RANGE];
  let mut left = [0; 2 * CODE_RANGE];
  let mut right = [0; 2 * CODE_RANGE];

  loop {
    let mut p0 = 0;
    let mut p1 = CODE_RANGE + 1;

    for i in 0..CODE_RANGE {
      tree_count[i] = entries[i].count;
    }

    for i in CODE_RANGE..(2 * CODE_RANGE) {
      tree_count[i] = 0xFFFFFFFF;
    }

    for w in (CODE_RANGE + 1)..(2 * CODE_RANGE) {
      if tree_count[p0] <= tree_count[p1] {
        left[w] = p0;
        p0 += 1;
      } else {
        left[w] = p1;
        p1 += 1;
      }

      if tree_count[p0] <= tree_count[p1] {
        right[w] = p0;
        p0 += 1;
      } else {
        right[w] = p1;
        p1 += 1;
      }

      tree_count[w] = tree_count[left[w]] + tree_count[right[w]];
    }

    bit_len[2 * CODE_RANGE - 1] = 0;
    for i in ((CODE_RANGE + 1)..(2 * CODE_RANGE)).rev() {
      bit_len[left[i]] = bit_len[i] + 1;
      bit_len[right[i]] = bit_len[i] + 1;
    }

    let mut max_len: u32 = 0;
    for i in 0..CODE_RANGE {
      max_len = u32::max(bit_len[i], max_len);
      entries[i].code_len = bit_len[i];
    }

    if max_len <= CODE_LENGTH_LIMIT as u32 {
      break;
    }

    let bs = max_len - (CODE_LENGTH_LIMIT as u32);
    for i in 0..CODE_RANGE {
      entries[i].count >>= bs;
      if entries[i].count == 0 {
        entries[i].count = 1;
      }
    }
  }

  let mut res: Vec<CodingTableEntry> = Vec::new();
  for _ in 0..CODE_RANGE {
    res.push(CodingTableEntry {
      code: 0,
      code_len: 0,
    });
  }

  for i in 0..CODE_RANGE {
    res[entries[i].sym as usize].code_len = entries[i].code_len;
  }

  assign_codes(&mut res);

  res
}

pub struct DecodingTable {
  pub base_code: Vec<u32>,
  pub base_offset: Vec<u32>,
  pub sym_table: Vec<u32>,
}

pub fn create_decode_table(frame: &[u8]) -> DecodingTable {
  let table = &create_code_table(frame);
  let (base_code, base_offset) = left_assign_codes(table, CODE_LENGTH_LIMIT);

  let mut sym_table = vec![0 as u32; CODE_RANGE];
  let mut cur_offset = base_offset.clone();

  for i in 0..CODE_RANGE {
    let len = table[i].code_len;
    sym_table[cur_offset[len as usize] as usize] = i as u32;
    cur_offset[len as usize] += 1;
  }

  DecodingTable {
    base_code,
    base_offset,
    sym_table,
  }
}
