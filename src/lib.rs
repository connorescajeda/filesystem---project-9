#![cfg_attr(not(test), no_std)]

#[derive(Copy, Clone, PartialEq, Eq)]
pub enum FileSystemResult<T: Copy + Clone> {
    Ok(T),
    Err(FileSystemError),
}

impl<T: Copy + Clone> FileSystemResult<T> {
    pub fn unwrap(&self) -> T {
        match self {
            FileSystemResult::Ok(v) => *v,
            FileSystemResult::Err(e) => panic!("Error: {e:?}"),
        }
    }
}

#[derive(Copy, Clone, Eq, PartialEq, Debug)]
pub enum FileSystemError {
    FileNotFound,
    FileNotOpen,
    NotOpenForRead,
    NotOpenForWrite,
    TooManyOpen,
    TooManyFiles,
    AlreadyOpen,
    DiskFull,
    FileTooBig,
    FilenameTooLong,
}

#[derive(Debug, Copy, Clone)]
pub struct FileInfo<const MAX_BLOCKS: usize, const BLOCK_SIZE: usize> {
    inode: Inode<MAX_BLOCKS, BLOCK_SIZE>,
    inode_num: usize,
    current_block: usize,
    offset: usize,
    writing: bool,
    reading: bool,
    block_buffer: [u8; BLOCK_SIZE],
}

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub struct Inode<const MAX_BLOCKS: usize, const BLOCK_SIZE: usize> {
    bytes_stored: u16,
    blocks: [u8; MAX_BLOCKS],
}

const INODE_FULL_BLOCK: usize = 0;
const DATA_FULL_BLOCK: usize = INODE_FULL_BLOCK + 1;
const INODE_TABLE_START: usize = DATA_FULL_BLOCK + 1;
//bitwise OR
//buffer[0] |= 1;
//buffer[0] |= 2 ;
//buffer[0] |= 4;
//buffer[0] |= (1 << bit)
// inode19
// 19/8 =2
//19 % 8 = 3
//buffer[2] = (1<< 3)
#[derive(core::fmt::Debug)]
pub struct FileSystem<
    const MAX_OPEN: usize,
    const BLOCK_SIZE: usize,
    const NUM_BLOCKS: usize,
    const MAX_FILE_BLOCKS: usize,
    const MAX_FILE_BYTES: usize,
    const MAX_FILES_STORED: usize,
    const MAX_FILENAME_BYTES: usize,
> {
    open: [Option<FileInfo<MAX_FILE_BLOCKS, BLOCK_SIZE>>; MAX_OPEN],
    disk: ramdisk::RamDisk<BLOCK_SIZE, NUM_BLOCKS>,
    block_buffer: [u8; BLOCK_SIZE],
    file_content_buffer: [u8; MAX_FILE_BYTES],
    directory_buffer: [u8; MAX_FILE_BYTES],
    open_inodes: [bool; MAX_FILES_STORED],
}

impl<
        const MAX_OPEN: usize,
        const BLOCK_SIZE: usize,
        const NUM_BLOCKS: usize,
        const MAX_FILE_BLOCKS: usize,
        const MAX_FILE_BYTES: usize,
        const MAX_FILES_STORED: usize,
        const MAX_FILENAME_BYTES: usize,
    >
    FileSystem<
        MAX_OPEN,
        BLOCK_SIZE,
        NUM_BLOCKS,
        MAX_FILE_BLOCKS,
        MAX_FILE_BYTES,
        MAX_FILES_STORED,
        MAX_FILENAME_BYTES,
    >
{
    pub fn new(disk: ramdisk::RamDisk<BLOCK_SIZE, NUM_BLOCKS>) -> Self {
        assert_eq!(MAX_FILE_BYTES, MAX_FILE_BLOCKS * BLOCK_SIZE);
        assert!(NUM_BLOCKS <= u8::MAX as usize);
        assert!(MAX_FILE_BYTES <= u16::MAX as usize);
        let block_bits = BLOCK_SIZE * 8;
        assert!(MAX_FILES_STORED <= block_bits);
        assert!(MAX_FILES_STORED <= u16::MAX as usize);
        let result = Self {
            open: [None; MAX_OPEN],
            disk,
            block_buffer: [0; BLOCK_SIZE],
            file_content_buffer: [0; MAX_FILE_BYTES],
            directory_buffer: [0; MAX_FILE_BYTES],
            open_inodes: [false; MAX_FILES_STORED],
        };
        assert!(result.num_inode_blocks() * 2 < NUM_BLOCKS);
        assert!(result.num_data_blocks() <= block_bits);
        assert_eq!(
            result.num_data_blocks() + result.num_inode_blocks() + 2,
            NUM_BLOCKS
        );
        assert!(result.num_inode_entries() <= u16::MAX as usize);
        assert!(result.num_inode_blocks() <= MAX_FILE_BLOCKS);
        result
    }

    pub fn list_directory(&mut self) -> FileSystemResult<(usize, [[u8; MAX_FILENAME_BYTES]; MAX_FILES_STORED])> {
        self.get_directory();
        let mut count = 0;
        let mut filenames = [['\0' as u8; MAX_FILENAME_BYTES]; MAX_FILES_STORED];
        for (i, c) in self.directory_buffer.iter().enumerate() {
            if i % MAX_FILENAME_BYTES == 0 && *c != 0 as u8 {
                count += 1;
                filenames[count - 1][i % MAX_FILENAME_BYTES] = *c;
            } else if i % MAX_FILENAME_BYTES != 0{
                filenames[count - 1][i % MAX_FILENAME_BYTES] = *c;
            } else{
                break;
            }
        }

        return FileSystemResult::Ok((count, filenames))
    }

    pub fn max_file_size(&self) -> usize {
        MAX_FILE_BLOCKS * BLOCK_SIZE
    }

    pub fn num_inode_bytes(&self) -> usize {
        2 + MAX_FILE_BLOCKS
    }

    pub fn inodes_per_block(&self) -> usize {
        BLOCK_SIZE / self.num_inode_bytes()
    }

    pub fn num_inode_blocks(&self) -> usize {
        MAX_FILES_STORED / self.inodes_per_block()
    }

    pub fn num_data_blocks(&self) -> usize {
        NUM_BLOCKS - self.num_inode_blocks() - 2
    }

    pub fn num_inode_entries(&self) -> usize {
        self.inodes_per_block() * self.num_inode_blocks() * self.num_inode_bytes()
    }

    pub fn first_data_block(&self) -> usize {
        2 + self.num_inode_blocks()
    }
    pub fn get_open(&self) -> [Option<FileInfo<MAX_FILE_BLOCKS, BLOCK_SIZE>>; MAX_OPEN]{
        self.open
    }
    pub fn get_file_content_buffer(&self) -> [u8; MAX_FILE_BYTES]{
        self.file_content_buffer
    }

    pub fn open_read(&mut self, filename: &str) -> FileSystemResult<usize> {
        self.get_directory();

        let mut namebuffer = ['\0'; MAX_FILENAME_BYTES];
        for (i, c) in filename.chars().enumerate() {
            namebuffer[i] = c
        }
        let mut found = false;
        let mut name_count = 0;
        let mut char_count = 0;
        let mut name_flag = true;
        let mut ignore = false;
        let mut inode_num = 0;

        for (i, value) in self.directory_buffer.iter().enumerate() {
            if ignore {
                if i % MAX_FILENAME_BYTES == 0 {
                    ignore = false;
                    name_flag = true;
                    char_count = 0;
                    name_count += 1;
                }
                    
            } 
            if !ignore {
                if *value as u8 as char != namebuffer[char_count as usize % MAX_FILENAME_BYTES] {
                    ignore = true;
                }
                char_count += 1;

                if char_count == MAX_FILENAME_BYTES && name_flag {
                    name_count += 1;
                    found = true;
                    break
                }
            }
        }
        if found {
            inode_num = name_count;
            if self.open_inodes[inode_num] {
                return FileSystemResult::Err(FileSystemError::AlreadyOpen);
            } else {
                let inode_start = inode_num * self.num_inode_bytes();
                let data = ((self.file_content_buffer[inode_start] as u16)<<8) | self.file_content_buffer[inode_start+1] as u16;
                let mut inode_blocks = [self.file_content_buffer[inode_start+2];MAX_FILE_BLOCKS];
                let mut count = 1;
                for i in inode_start+2..inode_start+self.num_inode_bytes(){
                    let block = self.file_content_buffer[i];
                    if inode_blocks.contains(&block) {

                    } else{
                        inode_blocks[count] = block;
                        count += 1;
                    }
                }
                let rebuilt_inode = Inode {
                    bytes_stored: data,
                    blocks: inode_blocks,
                };

                let mut buffer = [0; BLOCK_SIZE];
                self.disk.read(rebuilt_inode.blocks[0].into(), &mut buffer);

                //let mut file_offset = 0;
                
                // for (i, value) in buffer.iter().enumerate() {
                //     if *value == 0 as u8{
                //         file_offset = i;
                //         break;
                //     }
                // }
                
                let file_table_entry: FileInfo<MAX_FILE_BLOCKS, BLOCK_SIZE> = FileInfo{
                inode: rebuilt_inode,
                inode_num: inode_num,
                current_block: rebuilt_inode.blocks[0].into(),
                offset: 0,
                writing: false,
                reading: true,
                block_buffer: buffer
                };

                self.open_inodes[inode_num] = true;
                let mut fd = 0;
                for i in self.open{
                    if i.is_none() {
                        self.open[fd] = Some(file_table_entry);
                        break;
                    }
                    fd += 1
                }

                return FileSystemResult::Ok(fd)
            }
        } else {
            return FileSystemResult::Err(FileSystemError::FileNotFound);
        }
    }

    pub fn get_inode_table(&mut self){
        for i in 0..self.num_inode_blocks(){
            let mut buffer = [0; BLOCK_SIZE];
            self.disk.read(i+2, &mut buffer);
            for (j, value) in buffer.iter().enumerate() {
                self.file_content_buffer[j + (i * BLOCK_SIZE)] = *value;
            }
        }
    }

    pub fn get_directory(&mut self) {
        self.get_inode_table();
        let mut dir_blocks = [0; MAX_FILE_BLOCKS];
        let mut count = 0;
        for i in 2..self.num_inode_bytes() {
            let block = self.file_content_buffer[i];
            if dir_blocks.contains(&block){

            } else{
                dir_blocks[count] = block;
                count += 1;
            }
        }

        for (i, block) in dir_blocks.iter().enumerate() {
            if i == count {
                break;
            }
            let mut buffer = [0; BLOCK_SIZE];
            self.disk.read(*block as usize, &mut buffer);
            for (j, value) in buffer.iter().enumerate(){
                self.directory_buffer[j + (BLOCK_SIZE * i)] = *value;
            }
        }
        
    }

    pub fn open_create(&mut self, filename: &str) -> FileSystemResult<usize> {
        let mut itable = [0 ; BLOCK_SIZE];
        let mut datatable: [u8; BLOCK_SIZE] = [0 ; BLOCK_SIZE];
        let mut namebuffer = ['\0'; MAX_FILENAME_BYTES];
        if filename.len() > MAX_FILENAME_BYTES {
            return FileSystemResult::Err(FileSystemError::FilenameTooLong);
        }
        for (i, c) in filename.chars().enumerate() {
            namebuffer[i] = c
        }

        self.disk.read(INODE_FULL_BLOCK, &mut itable);
        self.disk.read(DATA_FULL_BLOCK, &mut datatable);

        if itable[0] & (1 << 0) == 0 {
            let num_active_blocks = 2 + self.num_inode_blocks();
            for i in 0..num_active_blocks + 1{
                let block = i / 8;
                let bit = i % 8;
                datatable[block] |= 1 << bit;
            }
            itable[0] = 1 << 0;
            self.disk.write(DATA_FULL_BLOCK, &mut datatable);
            self.disk.write(INODE_FULL_BLOCK, &mut itable);

            let data_block = 2 + self.num_inode_blocks();

            let inode:Inode<MAX_FILE_BLOCKS, BLOCK_SIZE> = Inode {
                bytes_stored : 0,
                blocks: [data_block.try_into().unwrap(); MAX_FILE_BLOCKS],
            };
            
            let piece1 = (inode.bytes_stored >> 8) as u8;
            let piece2 = inode.bytes_stored as u8;
            let mut ibuffer = [0;BLOCK_SIZE];
            self.disk.read(INODE_TABLE_START, &mut ibuffer);
            ibuffer[0] = piece1;
            ibuffer[1] = piece2;
            let mut start = 2;
            for block in inode.blocks{
                ibuffer[start] = block;
                start += 1;
            }
            self.disk.write(INODE_TABLE_START, &mut ibuffer);
        } 
        //Finding if file is new or has an inode
        let mut create_itable_buffer = [0 ; BLOCK_SIZE]; 
        let mut itable_count = 0;

        for i in 2..self.num_inode_blocks(){
            self.disk.read(i, &mut create_itable_buffer);
            for i in 0..BLOCK_SIZE {
                self.file_content_buffer[i + (itable_count * BLOCK_SIZE)] = create_itable_buffer[i]
            }
            itable_count += 1;
        }

        let block_count = MAX_FILE_BLOCKS;
        let mut directory_blocks = [0; MAX_FILE_BLOCKS];
        let mut count = 0;

        //Finds what blocks the directory is using
        for i in self.file_content_buffer {
            if count < block_count + 2{
                if count > 1 {
                    if !directory_blocks.contains(&i){
                        directory_blocks[count - 2] = i
                    }
                }
                count += 1
            } else{
                break;
            }
            
        }
        
        self.get_directory();

        let mut found = false;
        let mut name_count = 0;
        let mut char_count = 0;
        let mut name_flag = true;
        let mut ignore = false;
        let inode_num;

        for (i, value) in self.directory_buffer.iter().enumerate() {
            if ignore {
                if i % MAX_FILENAME_BYTES == 0 {
                    ignore = false;
                    name_flag = true;
                    char_count = 0;
                    name_count += 1;
                }
                    
            } 
            if !ignore {
                if *value as u8 as char != namebuffer[char_count as usize % MAX_FILENAME_BYTES] {
                    ignore = true;
                }
                char_count += 1;

                if char_count == MAX_FILENAME_BYTES && name_flag {
                    name_count += 1;
                    found = true;
                    break
                }
            }
        }

        if found {
            inode_num = name_count;       
        } else {
            inode_num = 0;
        }
        if inode_num != 0 {
            if self.open_inodes[inode_num] {
                return FileSystemResult::Err(FileSystemError::AlreadyOpen)
            } else {
                let inode_start = inode_num * self.num_inode_bytes();
                let mut data_table_buffer = [0; BLOCK_SIZE];
                self.disk.read(DATA_FULL_BLOCK, &mut data_table_buffer);
                let mut count = 0;
                let mut in_use = 0;
                for i in inode_start..inode_start + self.num_inode_bytes(){
                    if count < 2 {
                        self.file_content_buffer[i as usize] = 0 as u8;
                        } 
                    if count == 2{
                        in_use = self.file_content_buffer[i as usize];
                    }
                    if count > 2{
                        let block = i / 8;
                        let bit = i % 8;
                        data_table_buffer[block as usize] &= !(1 << bit);
                        self.file_content_buffer[i as usize] = in_use;
                    }
                    count += 1;
    
                }
                let mut start_block = 2;
    
                self.file_content_buffer = self.write_to_inode_table(start_block);
                
                self.disk.write(DATA_FULL_BLOCK, &mut data_table_buffer);
    
                let data = ((self.file_content_buffer[inode_start as usize] as u16)<<8) | self.file_content_buffer[inode_start as usize+1] as u16;
                let mut inode_blocks = [0;MAX_FILE_BLOCKS];
                    let mut count = 0;
                    for i in inode_start+2..inode_start+self.num_inode_bytes(){
                        if inode_blocks.contains(&(self.file_content_buffer[i as usize])) {
    
                        } else{
                            inode_blocks[count] = self.file_content_buffer[i as usize];
                            count += 1;
                        }
                    }
    
                let rebuilt_inode = Inode {
                    bytes_stored: data,
                    blocks: inode_blocks,
                };   
                let mut new_buffer = [0 as u8;BLOCK_SIZE]; 
                self.disk.read(rebuilt_inode.blocks[0].into(), &mut new_buffer);
    
                let file_table_entry: FileInfo<MAX_FILE_BLOCKS, BLOCK_SIZE> = FileInfo{
                    inode: rebuilt_inode,
                    inode_num: inode_num as usize,
                    current_block: rebuilt_inode.blocks[0].into(),
                    offset: 0,
                    writing: false,
                    reading: false,
                    block_buffer: new_buffer
                    };
    
                let mut fd = 0;
                for i in self.open{
                    if i.is_none(){
                        self.open[fd] = Some(file_table_entry);
                        break;
                    } 
                    fd+=1;
                    
                    
                }
                return FileSystemResult::Ok(fd);
            }
  
        
        } else {
            if datatable[BLOCK_SIZE - 1] == u8::MAX {
                return FileSystemResult::Err(FileSystemError::DiskFull);
            }
                self.disk.read(0, &mut itable);
                let inode_info = self.return_open_inode();
                if inode_info[2] == MAX_FILES_STORED as u8 {
                    return FileSystemResult::Err(FileSystemError::TooManyFiles)
                }
                let inode_num = inode_info[2];
                itable[inode_info[0] as usize] |= 1 << inode_info[1];
                self.disk.write(0, &mut itable);
                
                self.disk.read(1, &mut datatable);

                let mut blocks = [0; 2];
                let data_block = self.return_open_data();
                datatable[data_block[0] as usize] |= 1 << data_block[1];
                blocks[0] = data_block[2];

                let inode_start = inode_num as u16 * self.num_inode_bytes() as u16;
                let inode_table_block = 2 + (inode_start / BLOCK_SIZE as u16);

                // let data_block2;
                
                // if inode_start + self.num_inode_bytes() as u8 >= BLOCK_SIZE as u8 {
                //     block_flag = true;
                //     let data_block2 = self.return_open_data();
                //     datatable[data_block2[0] as usize] |= 1 << data_block2[1];
                //     blocks[1] = data_block2[2];
                // }

                self.disk.write(1, &mut datatable);

                let inode_blocks = [blocks[0]; MAX_FILE_BLOCKS];
                // if blocks[1] != 0 {
                //     inode_blocks[1] = blocks[1];
                // }

                let inode:Inode<MAX_FILE_BLOCKS, BLOCK_SIZE> = Inode {
                    bytes_stored: 0,
                    blocks: inode_blocks,
                };
                let piece1 = (inode.bytes_stored >> 8) as u8;
                let piece2 = inode.bytes_stored as u8;

                self.file_content_buffer[inode_start as usize] = piece1;
                self.file_content_buffer[inode_start as usize + 1] = piece2;

                let mut count = 2;
                for block in inode.blocks {
                    self.file_content_buffer[inode_start as usize + count as usize] = block;
                    count += 1;
                }
                self.file_content_buffer = self.write_to_inode_table( inode_table_block as usize);

                let directory_index = MAX_FILENAME_BYTES * (inode_num - 1) as usize;
                let index_start = directory_index % BLOCK_SIZE;
                
                let mut count = 0;
                
                if BLOCK_SIZE - index_start < MAX_FILENAME_BYTES {
                    let new_data_block = self.return_open_data();
                    datatable[new_data_block[0] as usize] |= 1 << new_data_block[1];
                    self.file_content_buffer = self.add_new_data_to_inode(0, new_data_block[2]);
                    self.disk.write(1, &datatable);
                }
                
                
                let mut dir_blocks = [0; MAX_FILE_BLOCKS];
                for i in 2..self.num_inode_bytes(){
                    let block = self.file_content_buffer[i];
                    if dir_blocks.contains(&block){

                    } else{
                        dir_blocks[count] = block;
                        count += 1;
                    }
                }
               
                
                for (i, block) in dir_blocks.iter().enumerate() {
                    if *block == 0 as u8 {
                        break;
                    }
                    let mut buffer = [0; BLOCK_SIZE];
                    self.disk.read(*block as usize, &mut buffer);
                    for (j, value) in buffer.iter().enumerate(){
                        self.directory_buffer[j + (BLOCK_SIZE * i)] = *value;
                    }
                }
                
                let mut _count = 0;
                for i in index_start..index_start + MAX_FILENAME_BYTES{
                    self.directory_buffer[i] = namebuffer[_count] as u8;
                    _count += 1;
                }
                
                //println!("{:?}", dir_blocks);
                self.directory_buffer = self.write_to_dir(dir_blocks);
                //println!("{:?}", self.directory_buffer);
                let mut block_buffer = [0 as u8;BLOCK_SIZE]; 
                self.disk.read(blocks[0] as usize, &mut block_buffer);

                let file_table_entry: FileInfo<MAX_FILE_BLOCKS, BLOCK_SIZE> = FileInfo{
                    inode: inode,
                    inode_num: inode_num as usize,
                    current_block: blocks[0].into(),
                    offset: 0,
                    writing: false,
                    reading: false,
                    block_buffer: block_buffer
                    };
                let mut fd = 0;
                for i in self.open{
                    if i.is_none(){
                        self.open[fd] = Some(file_table_entry);
                        break;
                    }
                    fd+=1;
                }
                self.open_inodes[inode_num as usize] = true;
                return FileSystemResult::Ok(fd);
            }

    }

    pub fn open_append(&mut self, filename: &str) -> FileSystemResult<usize> {
        self.get_directory();

        let mut namebuffer = ['\0'; MAX_FILENAME_BYTES];
        for (i, c) in filename.chars().enumerate() {
            namebuffer[i] = c 
        }
        let mut found = false;
        let mut name_count = 0;
        let mut char_count = 0;
        let mut name_flag = true;
        let mut ignore = false;
        let inode_num;

        for (i, value) in self.directory_buffer.iter().enumerate() {
            if ignore {
                if i % MAX_FILENAME_BYTES == 0 {
                    ignore = false;
                    name_flag = true;
                    char_count = 0;
                    name_count += 1;
                }
                    
            } 
            if !ignore {
                if *value as u8 as char != namebuffer[char_count as usize % MAX_FILENAME_BYTES] {
                    ignore = true;
                }
                char_count += 1;

                if char_count == MAX_FILENAME_BYTES && name_flag {
                    name_count += 1;
                    found = true;
                    break
                }
            }
        }
        if found {
            inode_num = name_count;
            if self.open_inodes[inode_num] {
                
                return FileSystemResult::Err(FileSystemError::AlreadyOpen);
            } else {
                let inode_start = inode_num * self.num_inode_bytes();
                let data = ((self.file_content_buffer[inode_start] as u16)<<8) | self.file_content_buffer[inode_start+1] as u16;
                let mut inode_blocks = [0;MAX_FILE_BLOCKS];
                let mut count = 0;
                for i in inode_start+2..inode_start+self.num_inode_bytes(){
                    if inode_blocks.contains(&self.file_content_buffer[i]) {

                    } else{
                        inode_blocks[count] = self.file_content_buffer[i];
                        count += 1;
                    }
                }
                let rebuilt_inode = Inode {
                    bytes_stored: data,
                    blocks: inode_blocks,
                };
                let mut buffer = [0; BLOCK_SIZE];
                let mut c = 0;
                for i in inode_blocks {
                    if i == 0 {
                        break;
                    } else {
                        c += 1;
                    }
                }
                let current: u8;
                if c != 0 {
                    self.disk.read(rebuilt_inode.blocks[c - 1].into(), &mut buffer);
                    current = rebuilt_inode.blocks[c - 1];
                } else {;
                    self.disk.read(rebuilt_inode.blocks[c].into(), &mut buffer);
                    current = rebuilt_inode.blocks[c].into();
                }
                

                let mut offset = 0;
                for (i, value) in buffer.iter().enumerate() {
                    if *value == (0 as u8){
                        offset = i;
                        break;
                    }
                }

                let file_table_entry: FileInfo<MAX_FILE_BLOCKS, BLOCK_SIZE> = FileInfo{
                inode: rebuilt_inode,
                inode_num: inode_num,
                current_block: current.into(),
                offset: offset,
                writing: true,
                reading: false,
                block_buffer: buffer
                };

                self.open_inodes[inode_num] = true;
                let mut fd = 0;
                for i in self.open{
                    if i.is_none() {
                        self.open[fd] = Some(file_table_entry);
                        break;
                    }
                    fd += 1;
                }
                return FileSystemResult::Ok(fd)
            }
        } else {
            return FileSystemResult::Err(FileSystemError::FileNotFound);
        }
    }
    
    pub fn read(&mut self, fd: usize, buffer: &mut [u8]) -> FileSystemResult<usize> {
        let mut og_file = self.open[fd];
        if og_file.is_none() {
            return FileSystemResult::Err(FileSystemError::FileNotOpen)    
        } else {
            let mut file = og_file.unwrap();
            if file.writing {
                return FileSystemResult::Err(FileSystemError::NotOpenForRead) 
            } else{
                file.reading = true;
                let start = file.inode_num * self.num_inode_bytes() + 2;
                let stop = file.inode_num * self.num_inode_bytes() + self.num_inode_bytes();
                let inode_blocks= &self.file_content_buffer[start..stop];
                let mut bytes: usize = 0;
                let mut block = file.current_block;
                let mut unique_blocks = [0 as usize; MAX_FILE_BLOCKS];
                let mut count = 1;
                let mut blocks_used = 0;
                
                for b in inode_blocks {
                    if unique_blocks.contains(&(*b as usize)) {

                    } else{
                        unique_blocks[count] = *b as usize;
                        if block == *b as usize {
                            blocks_used = count
                        }
                        count += 1;
                        if count == MAX_FILE_BLOCKS {
                            break;
                        }
                    }
                    
                }

                for j in 0..buffer.len() {     
                        
                    if blocks_used == count {
                        break;
                    }
                    let mut disk_buffer = [0; BLOCK_SIZE];
                    self.disk.read(block as usize, &mut disk_buffer);
                    if disk_buffer[file.offset] == 0{
                        break;
                    }
                    if file.offset == BLOCK_SIZE - 1{
                        if blocks_used == count {
                            break;
                        }  
                        buffer[j] = disk_buffer[file.offset];
                        blocks_used += 1;
                        
                        block = unique_blocks[blocks_used];
                        file.current_block = block as usize;
                        file.offset = 0;
                    } else {
                        buffer[j] = disk_buffer[file.offset];
                        file.offset += 1;
                    }
                    
                    bytes += 1;
                }

                self.open[fd] = Some(file);

            
                return FileSystemResult::Ok(bytes);
            }
        } 
    }
    
    pub fn write(&mut self, fd: usize, buffer: &[u8]) -> FileSystemResult<()> {
        let mut og_file = self.open[fd];
        let mut datatable = [0; BLOCK_SIZE];
        self.disk.read(1, &mut datatable);
        let max = 0;
        if og_file.is_none() {
            return FileSystemResult::Err(FileSystemError::FileNotOpen)    
        } else { 
            let mut file = og_file.unwrap(); 

            if file.reading {
                return FileSystemResult::Err(FileSystemError::NotOpenForWrite)  
            } else{
                let mut full = true;
                let num_blocks = NUM_BLOCKS / 8; // changed for magic number
                for (i, value) in datatable.iter().enumerate() {
                    //println!("{i}");
                    if i > num_blocks  {
                        break;
                    }
                    if *value != u8::MAX {
                        full = false;
                        break;
                    }
                }
                if full {
                    return FileSystemResult::Err(FileSystemError::DiskFull)
                }
                let mut full_flag = false;
                let mut blocks = [0;MAX_FILE_BLOCKS];
                let mut count = 0;
                let inode_start = file.inode_num * self.num_inode_bytes();
                for i in inode_start + 2.. inode_start + self.num_inode_bytes(){
                    if blocks.contains(&self.file_content_buffer[i]) {

                    } else{
                        if count >= MAX_FILE_BLOCKS {
                            return FileSystemResult::Err(FileSystemError::FileTooBig)
                        } 
                        blocks[count] = self.file_content_buffer[i];
                        count += 1;
                    }
                }
                if blocks.contains(&0) {

                } else{
                    full_flag = true;
                }

                file.writing = true;

                for (i, value) in buffer.iter().enumerate(){
                    if full_flag {
                        let mut buff = [0; BLOCK_SIZE];
                        self.disk.read(blocks[MAX_FILE_BLOCKS - 1] as usize, &mut buff);
                        if buff[BLOCK_SIZE - 2] != 0 {
                            
                            return FileSystemResult::Err(FileSystemError::FileTooBig)
                        }
                      
                    };
                    file.block_buffer[file.offset] = *value;
                    file.offset += 1;
                    if file.offset >= BLOCK_SIZE{
                        self.disk.write(file.current_block, &file.block_buffer);
                        let new_block = self.return_open_data();
                        datatable[new_block[0] as usize] |= 1 << new_block[1];
                        self.disk.write(1, &datatable);
                        self.file_content_buffer = self.add_new_data_to_inode(file.inode_num as u8, new_block[2]);
                        file.current_block = new_block[2] as usize;
                        file.offset = 0;
                        file.block_buffer = [0; BLOCK_SIZE];
                    }
                }
                og_file = Some(file);
                
                self.disk.write(file.current_block, &file.block_buffer);
        }
        
        self.open[fd] = og_file;
        return FileSystemResult::Ok(())
        }
            
    } 

    pub fn close(&mut self, fd: usize) -> FileSystemResult<()> {
        let file = self.open[fd];
        if file.is_none() {
            return FileSystemResult::Err(FileSystemError::FileNotOpen)    
        } else {
            let mut file = file.unwrap();
            self.get_inode_table();
            let inode_start = file.inode_num * self.num_inode_bytes();
            let mut used_blocks = [self.file_content_buffer[inode_start + 2] as u8;MAX_FILE_BLOCKS];
            let mut count = 1;
            for i in (inode_start + 2)..(inode_start + self.num_inode_bytes()){
                let value = self.file_content_buffer[i];
                if used_blocks.contains(&value) {
                } else{
                    used_blocks[count] = value;
                    count += 1;
                }
            }
            let mut total_bytes = 0;
           
            for i in 0..count {
                let mut buffer = [0; BLOCK_SIZE];
                self.disk.read(used_blocks[i] as usize, &mut buffer);
                
                for j in buffer {
                    if j != 0 {
                        total_bytes += 1;
                    } else{
                        break;
                    }
                }
            }
            
            let piece1 = (total_bytes >> 8) as u8;
            let piece2 = total_bytes as u8;
            self.file_content_buffer[inode_start] = piece1;
            self.file_content_buffer[inode_start + 1] = piece2;
            self.file_content_buffer = self.write_to_inode_table(2 + (inode_start / BLOCK_SIZE));
            self.open[fd] = None;
            self.open_inodes[file.inode_num] = false;
            return FileSystemResult::Ok(());
        }
        
    }
   
   pub fn add_new_data_to_inode(&mut self, inode_num: u8, new_data_block : u8) -> [u8;MAX_FILE_BYTES] {
        let mut buffer = [0; BLOCK_SIZE];
        let inode_start = inode_num as u16 * self.num_inode_bytes() as u16;
        let mut last_block = 0;
        let mut index = (inode_start + 3) as u16;
        let mut flag = false;
        for i in (inode_start + 2) as u16..(inode_start as u16 + self.num_inode_bytes() as u16){
            if !flag{
                flag = true;
                last_block = self.file_content_buffer[i as usize];
            } else if flag && self.file_content_buffer[i as usize] == last_block{
                break;
            } else if flag {
                last_block = self.file_content_buffer[i as usize];
                index = i as u16;
            }
        }
        self.file_content_buffer[index as usize] = new_data_block;
        self.write_to_inode_table(2 + inode_start as usize / BLOCK_SIZE)
        
    }

    pub fn write_to_dir(&mut self, dir_blocks: [u8; MAX_FILE_BLOCKS]) -> [u8;MAX_FILE_BYTES]{
        let mut dir_table_buffer = [0; BLOCK_SIZE];
        
        let mut blocks = [dir_blocks[0]; MAX_FILE_BLOCKS];
        let mut num_of_blocks = 1;
        for i in dir_blocks {
            if blocks.contains(&i) {

            } else{
                blocks[num_of_blocks] = i;
                num_of_blocks += 1;
            }
        }
        let mut block = blocks[0];
        let mut blocks_used = 0;
        let mut count = 0;
        for i in self.directory_buffer {   
            if block == 0  {
                break;
            } else {
                if BLOCK_SIZE - 1 == count{
                    self.disk.write(blocks[blocks_used].into(), &mut dir_table_buffer);
                    dir_table_buffer = [0; BLOCK_SIZE];
                    count = 0;
                    blocks_used += 1;
                    block = blocks[blocks_used];
                } else {
                    dir_table_buffer[count] = i;
                }
                count += 1
            }
            
        }
        return self.directory_buffer
        
    }

    pub fn write_to_inode_table(&mut self, start_block: usize)  -> [u8;MAX_FILE_BYTES]{
        let mut count = 0;
        let mut inode_table_buffer = [0; BLOCK_SIZE];
        let mut start = start_block;
        let mut break_flag = false;
        let total = MAX_FILES_STORED * self.num_inode_bytes();
        let mut total_count = 0;
        for i in self.file_content_buffer { 
            if total_count >= total {
                break;
            }
            if count + 1 == BLOCK_SIZE {
                self.disk.write(start, &mut inode_table_buffer);
                inode_table_buffer = [0; BLOCK_SIZE];
                count = 0;
                start += 1;
            } else {
                inode_table_buffer[count] = i;
            }
            count += 1;
            total_count += 1;
        }
        return self.file_content_buffer;
    }

    pub fn return_open_inode(&self) -> [u8; 3] {
        //itable[0] & (1 << 0) == 0 
        let mut buffer = [0; BLOCK_SIZE];
        self.disk.read(0, &mut buffer);
        let mut count = 0;
        let mut block_bit:[u8; 3] = [0 as u8; 3];
        for i in 0..BLOCK_SIZE {
            for j in 0..8{
                if buffer[i] & (1 << j) == 0 {
                    block_bit[0] = i as u8;
                    block_bit[1] = j;
                    block_bit[2] = count;
                    return block_bit;
                }
                count += 1;
            }
        }
        //NEED TO HANDLE AN ERROR HERE
        return [0,0,0];
    }

    pub fn return_open_data(&self) -> [u8; 3] {
        //itable[0] & (1 << 0) == 0 
        let mut buffer = [0; BLOCK_SIZE];
        self.disk.read(1, &mut buffer);
        let mut count = 0;
        let mut block_bit:[u8; 3] = [0 as u8; 3];
        for i in 0..BLOCK_SIZE {
            for j in 0..8{
                if buffer[i] & (1 << j) == 0 {
                    block_bit[0] = i as u8;
                    block_bit[1] = j;
                    block_bit[2] = count;
                    return block_bit;
                }
                count += 1;
            }
        }
        //NEED TO HANDLE AN ERROR HERE
        return block_bit;
    }

    
}

    



#[cfg(test)]
mod tests {
    use core::num;

    use super::*;

    const BLOCK_SIZE: usize = 64;
    const MAX_FILES_STORED: usize = 32;

    fn make_small_fs() -> FileSystem<16, 64, 255, 8, 512, 32, 8> {
        FileSystem::new(ramdisk::RamDisk::new())
    }

    #[test]
    fn test_short_write() {
        let mut sys = make_small_fs();
        let f1 = sys.open_create("one.txt").unwrap();
        sys.write(f1, "This is a test.".as_bytes()).unwrap();
        let mut buffer = [0; 50];
        sys.close(f1).unwrap();
        let f2 = sys.open_read("one.txt").unwrap();
        let bytes_read = sys.read(f2, &mut buffer).unwrap();
        assert_eq!(bytes_read, 15);
        let s = core::str::from_utf8(&buffer[0..bytes_read]).unwrap();
        assert_eq!(s, "This is a test.");
    }
        
    const LONG_DATA: &str = "This is a much, much longer message.
    It crosses a number of different lines in the text editor, all synthesized
    with the goal of exceeding the 64 byte block limit by a considerable amount.
    To that end, this text contains considerable excessive verbiage.";

    #[test]
    fn test_long_write() {
        assert_eq!(265, LONG_DATA.len());
        let mut sys = make_small_fs();
        let f1 = sys.open_create("one.txt").unwrap();
       // println!("{:?}", sys.disk);
        sys.write(f1, LONG_DATA.as_bytes()).unwrap();
        //println!("{:?}", sys.disk);
        sys.close(f1);
        
        let read = read_to_string(&mut sys, "one.txt");
        assert_eq!(read.as_str(), LONG_DATA);
    }

    fn read_to_string(
        sys: &mut FileSystem<16, BLOCK_SIZE, 255, 8, 512, 32, 8>,
        filename: &str,
    ) -> String {
        let fd = sys.open_read(filename).unwrap();
        let mut read = String::new();
        let mut buffer = [0; 10];
        loop {
            let num_bytes = sys.read(fd, &mut buffer).unwrap();
            let s = core::str::from_utf8(&buffer[0..num_bytes]).unwrap();
            read.push_str(s);
            if num_bytes < buffer.len() {
                sys.close(fd).unwrap();
                return read;
            }
        }
    }

    #[test]
    fn test_complex_1() {
        let one = "This is a message, a short message, but an increasingly long message.
        This is a message, a short message, but an increasingly long message.";
        let two = "This is the second message I have chosen to undertake in this particular test.
        This is a continuation of this ever-so-controversial second message.\n";
        let mut sys = make_small_fs();
        let f1 = sys.open_create("one.txt").unwrap();
        sys.write(f1, one[0..one.len() / 2].as_bytes()).unwrap();
        let f2 = sys.open_create("two.txt").unwrap();
        sys.write(f2, two[0..two.len() / 2].as_bytes()).unwrap();
        sys.write(f1, one[one.len() / 2..one.len()].as_bytes())
            .unwrap();
        sys.write(f2, two[two.len() / 2..two.len()].as_bytes())
            .unwrap();
        sys.close(f1).unwrap();
        sys.close(f2).unwrap();
        assert_eq!(one, read_to_string(&mut sys, "one.txt").as_str());
        assert_eq!(two, read_to_string(&mut sys, "two.txt").as_str());
    }

    #[test]
    fn test_complex_2() {
        let one = "This is a message, a short message, but an increasingly long message.
        This is a message, a short message, but an increasingly long message.";
        let two = "This is the second message I have chosen to undertake in this particular test.
        This is a continuation of this ever-so-controversial second message.\n";
        let mut sys = make_small_fs();
        let f1 = sys.open_create("one.txt").unwrap();
        sys.write(f1, one[0..one.len() / 2].as_bytes()).unwrap();

        let f2 = sys.open_create("two.txt").unwrap();
        //
        sys.write(f2, two[0..two.len() / 2].as_bytes()).unwrap();
        sys.close(f1).unwrap();
  
        sys.close(f2).unwrap();
        let f3 = sys.open_append("two.txt").unwrap();

        let f4 = sys.open_append("one.txt").unwrap();
        println!("{:?}", sys.list_directory().unwrap());
        sys.write(f4, one[one.len() / 2..one.len()].as_bytes())
            .unwrap();
        sys.write(f3, two[two.len() / 2..two.len()].as_bytes())
            .unwrap();
        
        sys.close(f1).unwrap();
        sys.close(f2).unwrap();
        
        assert_eq!(one, read_to_string(&mut sys, "one.txt").as_str());
        assert_eq!(two, read_to_string(&mut sys, "two.txt").as_str());
    }
    
    #[test]
    fn test_complex_3() {
        let one = "This is a message, a short message, but an increasingly long message.
        This is a message, a short message, but an increasingly long message.";
        let two = "This is the second message I have chosen to undertake in this particular test.
        This is a continuation of this ever-so-controversial second message.\n";
        let mut sys = make_small_fs();
        let f1 = sys.open_create("one.txt").unwrap();
        sys.write(f1, one.as_bytes()).unwrap();
        sys.close(f1).unwrap();
        println!("{:?}", sys.list_directory().unwrap().1);
        let f2 = sys.open_create("one.txt").unwrap();
        sys.write(f2, two.as_bytes()).unwrap();
        sys.close(f2).unwrap();

        assert_eq!(two, read_to_string(&mut sys, "one.txt").as_str());
    }

    #[test]
    fn test_file_not_found() {
        let mut sys = make_small_fs();
        let f1 = sys.open_create("one.txt").unwrap();
        sys.write(f1, "This is a test.".as_bytes()).unwrap();
        sys.close(f1).unwrap();
        match sys.open_read("one.tx") {
            FileSystemResult::Ok(_) => panic!("Shouldn't have found the file"),
            FileSystemResult::Err(e) => assert_eq!(e, FileSystemError::FileNotFound),
        }
    }

    #[test]
    fn test_file_not_open() {
        let mut sys = make_small_fs();
        let f1 = sys.open_create("one.txt").unwrap();
        sys.write(f1, "This is a test.".as_bytes()).unwrap();
        sys.close(f1).unwrap();
        let fd = sys.open_read("one.txt").unwrap();
        let mut buffer = [0; 10];
        match sys.read(fd + 1, &mut buffer) {
            FileSystemResult::Ok(_) => panic!("Should be an error!"),
            FileSystemResult::Err(e) => assert_eq!(e, FileSystemError::FileNotOpen),
        }
    }

    #[test]
    fn test_not_open_for_read() {
        let mut sys = make_small_fs();
        let f1 = sys.open_create("one.txt").unwrap();
        sys.write(f1, "This is a test.".as_bytes()).unwrap();
        let mut buffer = [0; 10];
        match sys.read(f1, &mut buffer) {
            FileSystemResult::Ok(_) => panic!("Should not work!"),
            FileSystemResult::Err(e) => assert_eq!(e, FileSystemError::NotOpenForRead),
        }
    }

    #[test]
    fn test_not_open_for_write() {
        let mut sys = make_small_fs();
        let f1 = sys.open_create("one.txt").unwrap();
        sys.write(f1, "This is a test.".as_bytes()).unwrap();
        sys.close(f1).unwrap();
        let f2 = sys.open_read("one.txt").unwrap();
        match sys.write(f2, "this is also a test".as_bytes()) {
            FileSystemResult::Ok(_) => panic!("Should be an error"),
            FileSystemResult::Err(e) => assert_eq!(e, FileSystemError::NotOpenForWrite),
        }
    }

    #[test]
    fn test_filename_too_long() {
        let mut sys = make_small_fs();
        match sys.open_create("this_is_an_exceedingly_long_filename_to_use.txt") {
            FileSystemResult::Ok(_) => panic!("This should be an error"),
            FileSystemResult::Err(e) => assert_eq!(e, FileSystemError::FilenameTooLong),
        }
    }

    #[test]
    fn test_already_open() {
        let mut sys = make_small_fs();
        let f1 = sys.open_create("one.txt").unwrap();
        sys.write(f1, "This is a test.".as_bytes()).unwrap();
        match sys.open_read("one.txt") {
            FileSystemResult::Ok(_) => panic!("Should be an error"),
            FileSystemResult::Err(e) => assert_eq!(e, FileSystemError::AlreadyOpen),
        }
    }

    #[test]
    fn test_file_too_big() {
        let mut sys = make_small_fs();
        let f1 = sys.open_create("one.txt").unwrap();
        for _ in 0..sys.max_file_size() - 1 {
            sys.write(f1, "A".as_bytes()).unwrap();
        }
        match sys.write(f1, "B".as_bytes()) {
            FileSystemResult::Ok(_) => panic!("Should be an error!"),
            FileSystemResult::Err(e) => assert_eq!(e, FileSystemError::FileTooBig),
        }
    }

    #[test]
    fn test_too_many_files() {
        let mut sys = make_small_fs();
        for i in 0..MAX_FILES_STORED - 1 {
            let filename = format!("file{i}");
            let f = sys.open_create(filename.as_str()).unwrap();
            let content = format!("This is sentence {i}");
            sys.write(f, content.as_bytes()).unwrap();
            sys.close(f).unwrap();
        }
        
        match sys.open_create("Final") {
            FileSystemResult::Ok(_) => panic!("This should be an error!"),
            FileSystemResult::Err(e) => assert_eq!(e, FileSystemError::TooManyFiles),
        }
    }

    #[test]
    fn test_disk_full() {
        let mut sys = make_small_fs();
        for i in 0..MAX_FILES_STORED - 1 {
            let filename = format!("file{i}");
            let f = sys.open_create(filename.as_str()).unwrap();
            for j in 0..sys.max_file_size() - 1 {
                match sys.write(f, "A".as_bytes()) {
                    FileSystemResult::Ok(_) => {}
                    FileSystemResult::Err(e) => {
                        assert_eq!(i, 30);
                        assert_eq!(j, 191);
                        assert_eq!(e, FileSystemError::DiskFull);
                        return;
                    }
                }
            }
            sys.close(f).unwrap();
        }
        panic!("The disk should have been full!");
    }
}