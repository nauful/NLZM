This project implements the paper titled "On the Implementation of Minimum Redundancy Prefix Codes" by Moffat and Turpin. The main idea is to switch representations to a specialization of arithmetic ranges rather than a tree. This allows constant-time decompression without constructing any lookup tables. The encoder is also adaptive rather than static, so code lengths do not need to be transmitted.