/*
 * This code is provided solely for the personal and private use of students
 * taking the CSC369H course at the University of Toronto. Copying for purposes
 * other than this use is expressly prohibited. All forms of distribution of
 * this code, including but not limited to public repositories on GitHub,
 * GitLab, Bitbucket, or any other online platform, whether as given or with
 * any changes, are expressly prohibited.
 *
 * Authors: Alexey Khrabrov, Karen Reid
 *
 * All of the files in this directory and all subdirectories are:
 * Copyright (c) 2019 Karen Reid
 */

/**
 * CSC369 Assignment 1 - a1fs driver implementation.
 */

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>

// Using 2.9.x FUSE API
#define FUSE_USE_VERSION 29
#include <fuse.h>

#include "a1fs.h"
#include "fs_ctx.h"
#include "options.h"
#include "map.h"

#define MAX(x, y) (((x) > (y)) ? (x) : (y))

//NOTE: All path arguments are absolute paths within the a1fs file system and
// start with a '/' that corresponds to the a1fs root directory.
//
// For example, if a1fs is mounted at "~/my_csc369_repo/a1b/mnt/", the path to a
// file at "~/my_csc369_repo/a1b/mnt/dir/file" (as seen by the OS) will be
// passed to FUSE callbacks as "/dir/file".
//
// Paths to directories (except for the root directory - "/") do not end in a
// trailing '/'. For example, "~/my_csc369_repo/a1b/mnt/dir/" will be passed to
// FUSE callbacks as "/dir".


/**
 * Initialize the file system.
 *
 * Called when the file system is mounted. NOTE: we are not using the FUSE
 * init() callback since it doesn't support returning errors. This function must
 * be called explicitly before fuse_main().
 *
 * @param fs    file system context to initialize.
 * @param opts  command line options.
 * @return      true on success; false on failure.
 */
static bool a1fs_init(fs_ctx *fs, a1fs_opts *opts)
{
	// Nothing to initialize if only printing help
	if (opts->help) return true;

	size_t size;
	void *image = map_file(opts->img_path, A1FS_BLOCK_SIZE, &size);
	if (!image) return false;

	return fs_ctx_init(fs, image, size);
}

/**
 * Cleanup the file system.
 *
 * Called when the file system is unmounted. Must cleanup all the resources
 * created in a1fs_init().
 */
static void a1fs_destroy(void *ctx)
{
	fs_ctx *fs = (fs_ctx*)ctx;
	if (fs->image) {
		munmap(fs->image, fs->size);
		fs_ctx_destroy(fs);
	}
}

/** Get file system context. */
static fs_ctx *get_fs(void)
{
	return (fs_ctx*)fuse_get_context()->private_data;
}

//HELPER FUNCTIONS -------------------------------------------------------------------------------------------------------------------------
//A function to flip a bit in a bitmap. K is the kth bit to flip.
void  flipBit(a1fs_blk_t* bitmap, int k)
{
	int i = k / 32;        //gives the corresponding index in the array A
	int pos = k % 32;      //gives the corresponding bit position in A[i]

	unsigned int flag = 1;   // flag = 0000.....00001
	flag = flag << pos;      // flag = 0000...010...000   (shifted k positions)

	bitmap[i] = bitmap[i] ^ flag;      // Flip the bit at the k-th position in A[i]
}


//Finds the first free entry in a block or inode bitmap. It returns the index of the first bit that is 0 (-1 on failure). totalBlocks is superblock->blocks_count
//Later to find the block that the index returned refers to you must do superblock->data_blocks + findZeroBit().
int findZeroBit(a1fs_blk_t* bitmap, a1fs_blk_t totalBlocks) {
	for (a1fs_blk_t i = 0; i < totalBlocks/32 + 1; i++) {
		if (bitmap[i] != 2147483647) {
			for (a1fs_blk_t k = 0; k < 32; k++) {//gives the corresponding bit position in A[i]
				unsigned int flag = 1;   // flag = 0000.....00001
				flag = flag << k;      // flag = 0000...010...000   (shifted k positions)
				if (i * 32 + k == totalBlocks) {
					return -1;
				}
				if ((bitmap[i] & flag) == 0) { //Then the bit at pos k is free
					a1fs_blk_t free = (i * 32 + k);
					return (int)(i * 32 + k);
				}
			}
		}
	}
	return -1;
}

//Checks if the kth bit is 0. Returns 0 on success, -1 on failure.
int checkZeroBit(a1fs_blk_t* bitmap, int k, a1fs_blk_t totalBlocks) {
	int i = k / 32;        //gives the corresponding index in the array A
	int pos = k % 32;      //gives the corresponding bit position in A[i]

	if (i * 32 + k >= (int)totalBlocks) { //Out of bounds
		return -1;
	}

	unsigned int flag = 1;   // flag = 0000.....00001
	flag = flag << pos;      // flag = 0000...010...000   (shifted k positions)

	if ((bitmap[i] & flag) == flag) {
		return -1;
	}
	return 0;
}

//Returns the index of the last valid extent in the table, size is the inode size.
int lastExtentIndex(int size, a1fs_extent_table* table) {
	int last_extent = 0; //The extent which contains the last block
	int blocks = (int)(size / A1FS_BLOCK_SIZE) + 1; //How many blocks this inode has allocated
	for (int i = 0; i < 512; i++) {
		blocks -= table->i_extents[i].count;
		if (blocks <= 0) {
			last_extent = i;
			break;
		}
	}
	return last_extent;
}

//Returns the index of the longest extend
int* longest_extent(a1fs_blk_t* bitmap, a1fs_blk_t totalBlocks) {
    int* toReturn = malloc(sizeof(int) * 2);
    toReturn[0] = -1;
    toReturn[1] = -1;
    int longest = 0;
    int best_start = 0;
    int current_length = 0;
    int current_best_start = 0;
    int is_on_extent = 0; //0 if the algorithm is not on an extent, 1 when algorithm is on an extent
    for (a1fs_blk_t i = 0; i < totalBlocks; i++) {
        if (bitmap[i] != 2147483647) {
            for (a1fs_blk_t k = 0; k < 32; k++) {//gives the corresponding bit position in A[i]
                unsigned int flag = 1;   // flag = 0000.....00001
                flag = flag << k;      // flag = 0000...010...000   (shifted k positions)
                if ((bitmap[i] & flag) == 0) { //Then the bit at pos k is free
                    if(is_on_extent == 0) { // if the algorithm just detected a free extent
                        is_on_extent = 1;
                        a1fs_blk_t free = (i * 32 + k);
                        current_best_start = free; // sets the start pos as the first point of the free extent
                        current_length = 1; // sets the length of the current free extent to 1
                    }
                    else{ //if the algorithm is already on a free extent
                        current_length += 1; //increase current extent length by one
                    }
                    if(current_length>longest){ //if the current extent length is bigger than the max length
                        longest = current_length; //set longest to the current extent length
                    }
                }
                else{
                    if(is_on_extent == 1) { // if the algorithm just got out of a free extent
                        if(longest == current_length){ // if this was the longest extent
                            best_start = current_best_start; // set best start to the current best start
                        }
                        current_length = 0; // set current length to 0
                    }
                    is_on_extent = 0; // since we are not on a extent anymore change this to 0
                }
            }
        }
    }
    if(longest > 0){
        toReturn[0] = best_start;
        toReturn[1] = longest;
        return toReturn;
    }
    return toReturn;
}

//A function to look up a path. It will return an array of ints where index 0 = the error type (0 on success, -1 for ENOTDIR, -2 for ENOENT), index 1 for the last existing inode number in path,
// index 2 the number of calls to strtok. Caller must deallocate, index 3 is the inode number of the parent of the last existing inode, index 4 denotes the index of the last inode in its parent
int* path_lookup(const char* path) {
	int* toReturn = malloc(sizeof(int) * 3);
	toReturn[0] = 0;
	toReturn[1] = 0;
	toReturn[2] = 0;
	toReturn[3] = 0;
	toReturn[4] = 0;
	// Get the array of inodes (inode table)
	fs_ctx * fs = get_fs();
	a1fs_superblock* superblock = (a1fs_superblock*)(fs->image + A1FS_BLOCK_SIZE);
	a1fs_blk_t inode_table_block = superblock->inode_table;
	a1fs_inode* inodes = (a1fs_inode*)(fs->image + (A1FS_BLOCK_SIZE * inode_table_block));
	int inode_number = 0;
	int parent_number = 0;


	//Get first token in path 
	char* path_copy = malloc((strlen(path)) * sizeof(char));
	strcpy(path_copy, path);
	char* dir_entry_name;
	const char delimiter[2] = "/";
	dir_entry_name = strtok(path_copy, delimiter); // Get first element in path 

	//Go to last token in path  
	while (dir_entry_name != NULL) {
		//Check if inode is a valid dir (if not switch index 0 of toReturn to -1 for ENOTDIR)
		if (((&inodes[inode_number])->mode & S_IFDIR) != S_IFDIR) {
			toReturn[0] = -1;
			toReturn[1] = inode_number;
			toReturn[3] = parent_number;
			free(path_copy);
			return toReturn;
		}
		//Some vars we need to keep track of
		int dir_entry_count = inodes[inode_number].size / sizeof(a1fs_dentry); // How many dir entries does the current inode have
		a1fs_extent_table* extent_table = (a1fs_extent_table*)(fs->image + (inodes[inode_number].extent_table * A1FS_BLOCK_SIZE)); //The array of extents corresponding to the inode
		int extent_number = 0;
		a1fs_dentry* dentry = (a1fs_dentry*)(fs->image + (A1FS_BLOCK_SIZE * extent_table->i_extents[0].start));
		//Loop over all dir entries
		int error = -2;  //Want to return this for ENOENT 
		int k = 0; //The index into a dir_entry block (this reflects changes when going to a new extent)
		for (int i = 0; i < dir_entry_count; i++) {
			if (strcmp(dentry[k].name, dir_entry_name) == 0) { //Is the dentry.name equal to the current token?
				toReturn[4] = i;
				parent_number = inode_number;
				inode_number = (int)dentry[k].ino;
				extent_number = 0;
				error = 0;
				break;
			}
			k++;
			if (k *sizeof(a1fs_dentry) == extent_table->i_extents[extent_number].count * A1FS_BLOCK_SIZE) {  //Is the next dir entry in a new extent?
				extent_number++;
				dentry = (a1fs_dentry*)(fs->image + (A1FS_BLOCK_SIZE * extent_table->i_extents[extent_number].start));
				k = 0;
				}
			}
		if (error != 0) {    
			toReturn[0] = -2;
			toReturn[1] = inode_number;
			toReturn[3] = parent_number;
			free(path_copy);
			return toReturn;
			}
		dir_entry_name = strtok(NULL, delimiter);
		toReturn[2] = toReturn[2] + 1;
	}
	free(path_copy);
	toReturn[1] = inode_number;
	toReturn[3] = parent_number;
	return toReturn;
}

char* create_string(const char *buf, off_t offset, int file_end_point){
    if(file_end_point >= (int) offset){
        return buf;
    }
    else{
        int zeros = (int) offset - file_end_point;
        int total_length = zeros + strlen(buf);
        char *new_buf = malloc(sizeof(char)*(total_length));
        memset(new_buf, 0 ,zeros);
        memcpy(&new_buf[zeros], buf, total_length);
        return new_buf;
    }
}

void free_extends(a1fs_inode *file, int size){
    fs_ctx *fs = get_fs();
    a1fs_superblock* superblock = (a1fs_superblock*)(fs->image + A1FS_BLOCK_SIZE);
    a1fs_blk_t* blockbitmap = (a1fs_blk_t*)(fs->image + (superblock->block_bitmap * A1FS_BLOCK_SIZE));
    a1fs_extent_table* file_extent_table = (a1fs_extent_table*)(fs->image + (file->extent_table * A1FS_BLOCK_SIZE));

    int toRemove = (file->size - size)/A1FS_BLOCK_SIZE; //The number of blocks to fully remove
    int extent_index = lastExtentIndex(file->size, file_extent_table); //Our last extent
    while (toRemove != 0) {
        int start = file_extent_table->i_extents[extent_index].start;
        int length = file_extent_table->i_extents[extent_index].count;
        int new_length = MAX(length - toRemove, 0); //How much of the current extent is left after trying to remove the remaining blocks to remove
        for (int i = length -1; i >= new_length; i--) { //We dealloc the blocks
            flipBit(blockbitmap, start + i);
        }
        toRemove -= (length - new_length);
        file_extent_table->i_extents[extent_index].count = new_length; //The extent has a new length
    }
    file->size -= size;
    if (file->size == 0) { //We have to dealloc the extent table
        flipBit(blockbitmap, file->extent_table);
    }
}


//HELPER END ---------------------------------------------------------------------------------------------------------------------


/**
 * Get file system statistics.
 *
 * Implements the statvfs() system call. See "man 2 statvfs" for details.
 * The f_bfree and f_bavail fields should be set to the same value.
 * The f_ffree and f_favail fields should be set to the same value.
 * The following fields can be ignored: f_fsid, f_flag.
 * All remaining fields are required.
 *
 * Errors: none
 *
 * @param path  path to any file in the file system. Can be ignored.
 * @param st    pointer to the struct statvfs that receives the result.
 * @return      0 on success; -errno on error.
 */
static int a1fs_statfs(const char *path, struct statvfs *st)
{
	(void)path;// unused
	fs_ctx *fs = get_fs();

	memset(st, 0, sizeof(*st));
	st->f_bsize   = A1FS_BLOCK_SIZE;
	st->f_frsize  = A1FS_BLOCK_SIZE;
	//TODO: fill in the rest of required fields based on the information stored
	// in the superblock
	a1fs_superblock* superblock = (a1fs_superblock*)(fs->image + A1FS_BLOCK_SIZE);
	st->f_blocks = superblock->blocks_count; /*size of fs in block units*/
	st->f_bfree = superblock->free_blocks_count;    /* # free blocks */
	st->f_bavail = superblock->free_blocks_count;   /* # free blocks for unprivileged users */
	st->f_files = superblock->inodes_count;    /* # inodes */
	st->f_ffree = superblock->free_inodes_count;    /* # free inodes */
	st->f_favail = superblock->free_inodes_count;;   /* # free inodes for unprivileged users */
	st->f_namemax = A1FS_NAME_MAX;

	return 0;
}

/**
 * Get file or directory attributes.
 *
 * Implements the lstat() system call. See "man 2 lstat" for details.
 * The following fields can be ignored: st_dev, st_ino, st_uid, st_gid, st_rdev,
 *                                      st_blksize, st_atim, st_ctim.
 * All remaining fields are required.
 *
 * NOTE: the st_blocks field is measured in 512-byte units (disk sectors);
 *       it should include any metadata blocks that are allocated to the 
 *       inode.
 *
 * NOTE2: the st_mode field must be set correctly for files and directories.
 *
 * Errors:
 *   ENAMETOOLONG  the path or one of its components is too long.
 *   ENOENT        a component of the path does not exist.
 *   ENOTDIR       a component of the path prefix is not a directory.
 *
 * @param path  path to a file or directory.
 * @param st    pointer to the struct stat that receives the result.
 * @return      0 on success; -errno on error;
 */
static int a1fs_getattr(const char *path, struct stat *st)
{
	if (strlen(path) >= A1FS_PATH_MAX) return -ENAMETOOLONG;
	memset(st, 0, sizeof(*st));
	//NOTE: This is just a placeholder that allows the file system to be mounted
	// without errors. You should remove this from your implementation.
	fs_ctx* fs = get_fs();
	a1fs_superblock* superblock = (a1fs_superblock*)(fs->image + A1FS_BLOCK_SIZE);
	a1fs_blk_t inode_table_block = superblock->inode_table;
	a1fs_inode* inodes = (a1fs_inode *)(fs->image + (A1FS_BLOCK_SIZE * inode_table_block));
	int* path_info = path_lookup(path);
	if (path_info[0] == -1) {
		return -ENOTDIR;
	}
	if (path_info[0] == -2) {
		return -ENOENT;
	}
	st->st_mode = inodes[path_info[1]].mode;
	st->st_nlink = inodes[path_info[1]].links;
	st->st_size = inodes[path_info[1]].size;
	st->st_blocks = inodes[path_info[1]].size / 512;
	st->st_mtim = inodes[path_info[1]].mtime;
	free(path_info);
	return 0;
}

/**
 * Read a directory.
 *
 * Implements the readdir() system call. Should call filler(buf, name, NULL, 0)
 * for each directory entry. See fuse.h in libfuse source code for details.
 *
 * Assumptions (already verified by FUSE using getattr() calls):
 *   "path" exists and is a directory.
 *
 * Errors:
 *   ENOMEM  not enough memory (e.g. a filler() call failed).
 *
 * @param path    path to the directory.
 * @param buf     buffer that receives the result.
 * @param filler  function that needs to be called for each directory entry.
 *                Pass 0 as offset (4th argument). 3rd argument can be NULL.
 * @param offset  unused.
 * @param fi      unused.
 * @return        0 on success; -errno on error.
 */
static int a1fs_readdir(const char *path, void *buf, fuse_fill_dir_t filler,
                        off_t offset, struct fuse_file_info *fi)
{
	(void)offset;// unused
	(void)fi;// unused
	fs_ctx *fs = get_fs();
	int* path_info = path_lookup(path);
	int inode_number = path_info[1];
	a1fs_superblock* superblock = (a1fs_superblock*)(fs->image + A1FS_BLOCK_SIZE);
	a1fs_blk_t inode_table_block = superblock->inode_table;
	a1fs_inode* inodes = (a1fs_inode*)(fs->image + (A1FS_BLOCK_SIZE * inode_table_block));

	int dir_entry_count = inodes[inode_number].size / sizeof(a1fs_dentry); // How many dir entries does the current inode have
	a1fs_extent_table* extent_table = (a1fs_extent_table*)(fs->image + (inodes[inode_number].extent_table * A1FS_BLOCK_SIZE)); //The array of extents corresponding to the inode
	int extent_number = 0;
	int k = 0; //Index inside of extent
	a1fs_dentry* dentry = (a1fs_dentry*)(fs->image + (A1FS_BLOCK_SIZE * extent_table->i_extents[0].start));

	//Do we still need this?
	if (filler(buf, ".", NULL, 0) != 0) { return -ENOMEM; }
	if (filler(buf, "..", NULL, 0) != 0) { return -ENOMEM; }

	for (int i = 0; i < dir_entry_count; i++) {
		if(filler(buf, dentry[k].name, NULL, 0) != 0) { return -ENOMEM; };
		k++;
		if (k * sizeof(a1fs_dentry) == extent_table->i_extents[extent_number].count * A1FS_BLOCK_SIZE) {  //Is the next dir entry in a new extent?
			extent_number++;
			dentry = (a1fs_dentry*)(fs->image + (A1FS_BLOCK_SIZE * extent_table->i_extents[extent_number].start));
			k = 0;
		}
	}
	free(path_info);
	return 0;
}


/**
 * Create a directory.
 *
 * Implements the mkdir() system call.
 *
 * Assumptions (already verified by FUSE using getattr() calls):
 *   "path" doesn't exist.
 *   The parent directory of "path" exists and is a directory.
 *   "path" and its components are not too long.
 *
 * Errors:
 *   ENOMEM  not enough memory (e.g. a malloc() call failed).
 *   ENOSPC  not enough free space in the file system.
 *
 * @param path  path to the directory to create.
 * @param mode  file mode bits.
 * @return      0 on success; -errno on error.
 */
static int a1fs_mkdir(const char* path, mode_t mode)
{
	mode = mode | S_IFDIR;
	fs_ctx* fs = get_fs();
	//TODO: create a directory at given path with given mode
	(void)path;
	(void)mode;
	a1fs_superblock* superblock = (a1fs_superblock*)(fs->image + A1FS_BLOCK_SIZE);
	a1fs_blk_t inode_table_block = superblock->inode_table;
	a1fs_inode* inodes = (a1fs_inode*)(fs->image + (A1FS_BLOCK_SIZE * inode_table_block));
	a1fs_blk_t* blockbitmap = (a1fs_blk_t*)(fs->image + (superblock->block_bitmap * A1FS_BLOCK_SIZE));
	a1fs_blk_t* inodebitmap = (a1fs_blk_t*)(fs->image + (superblock->inode_bitmap * A1FS_BLOCK_SIZE));
	a1fs_blk_t blocks_needed = 0;

	//First we have to check if there are enoguh free inodes and blocks. 
	int* path_info = path_lookup(path);
	//We begin by checking if the last existing inode (to which we add more dirs) has enough space for another a1fs_dir_entry
	if ((inodes[path_info[1]].size == 0)) {
		blocks_needed++;
	}
	if ((inodes[path_info[1]].size % A1FS_BLOCK_SIZE) < sizeof(a1fs_dentry)) {
		blocks_needed++;
	}

	//Now we reconstruct the path at which path_looup exited (that is the rest of the directories we want to create)
	char* path_copy = malloc((strlen(path)) * sizeof(char));
	if (path_copy == NULL) { return -ENOMEM; }
	strcpy(path_copy, path);
	const char delimiter[2] = "/";
	char* path_left;
	path_left = strtok(path_copy, delimiter);
	for (int i = 0; i < path_info[2]; i++) { path_left = strtok(NULL, delimiter); }


	if (superblock->free_blocks_count < blocks_needed || superblock->free_inodes_count < 1) {
		free(path_copy);
		return -ENOSPC;
	}
	//------------------------------------------------------------------------------------------------------
	//Now we have the rest of a path that we need to create dirs for. 
	char* dir_to_create = strtok(path_left, delimiter);

	//Parent dir vars
	a1fs_ino_t parent_inode_number = (a1fs_ino_t)path_info[1];
	a1fs_inode* parent = &inodes[parent_inode_number];

	//Child dir vars
	int i = findZeroBit(inodebitmap, superblock->inodes_count); //is non negative thx to previous check
	a1fs_ino_t child_inode_number = (a1fs_ino_t)(i);
	a1fs_inode* child = &inodes[child_inode_number];

	//Three cases. Case 0: parent is empty and doesnt have an extent block yet. Case 1:  we need to create a new extent, Case 2 there is enough space in current extent
	//Case 0
	if (parent->size == 0) {
		int i = findZeroBit(blockbitmap, superblock->blocks_count); //is non negative thx to previous check
		parent->extent_table = i;
		a1fs_extent_table* parent_extent_table = (a1fs_extent_table*)(fs->image + (parent->extent_table * A1FS_BLOCK_SIZE));
		flipBit(blockbitmap, i);

		//We need to set the values of the extent
		int k = findZeroBit(blockbitmap, superblock->blocks_count);
		a1fs_blk_t k2 = (a1fs_blk_t)k;
		parent_extent_table->i_extents[0].start = k2;
		parent_extent_table->i_extents[0].count = 1;
		flipBit(blockbitmap, k2);
		superblock->free_blocks_count -= 2;

	}

	int space_left = (int)(A1FS_BLOCK_SIZE - (parent->size % A1FS_BLOCK_SIZE));
	int space_used = (int)A1FS_BLOCK_SIZE - space_left;
	int extent_offset = 0;
	if (space_left == A1FS_BLOCK_SIZE && parent->size != 0) { //Case 1
		a1fs_extent_table* parent_extent_table = (a1fs_extent_table*)(fs->image + (parent->extent_table * A1FS_BLOCK_SIZE));
		int last_extent = lastExtentIndex((int)(parent->size), parent_extent_table);
		int last_block = parent_extent_table->i_extents[last_extent].start + parent_extent_table->i_extents[last_extent].count - 1;
		//We try to extend last last extent by one to keep fragmentation low
		int k = checkZeroBit(blockbitmap, last_block + 1, superblock->blocks_count);
		if (k == 0) {
			parent_extent_table->i_extents[last_extent].count += 1;
			flipBit(blockbitmap, last_block + 1);
		}
		else { // Need to create new extent
			if (last_extent + 1 >= 512) {
				free(path_copy);
				return -ENOSPC;
			}
			int i = findZeroBit(blockbitmap, superblock->blocks_count);
			printf("find zero bit in case 1: %d \n", i);
			parent_extent_table->i_extents[last_extent + 1].start = i;
			printf("extent start in case 1: %d \n", parent_extent_table->i_extents[last_extent + 1].start);
			parent_extent_table->i_extents[last_extent + 1].count = 1;
			flipBit(blockbitmap, i);
			extent_offset += 1;
		}

	}
	//Go to last position in extent table (using size of inode)
	a1fs_extent_table* parent_extent_table = (a1fs_extent_table*)(fs->image + (parent->extent_table * A1FS_BLOCK_SIZE));
	printf("case 2 \n");
	int last_extent = lastExtentIndex((int)(parent->size), parent_extent_table) + extent_offset;
	printf("last extent: %d \n", last_extent);
	printf("extent offset: %d \n", extent_offset);
	int offset = (int)(space_used);
	int last_block = (int)(parent_extent_table->i_extents[last_extent].start + parent_extent_table->i_extents[last_extent].count) - 1;
	printf("last extent start: %d \n", parent_extent_table->i_extents[last_extent].start);
	printf("offset: %d \n", offset);
	printf("last block: %d \n", last_block);

	//Now we can add the dir entry 
	a1fs_dentry* dentry = fs->image + (last_block * A1FS_BLOCK_SIZE) + offset;
	strncpy(dentry->name, dir_to_create, A1FS_NAME_MAX);
	dentry->ino = child_inode_number;

	//Now we update the parent's info
	parent->size = parent->size + sizeof(a1fs_dentry);
	parent->links = parent->links + 1;
	clock_gettime(CLOCK_REALTIME, &parent->mtime);

	//Now we update the child's info
	child->size = 0;
	child->links = 2;
	clock_gettime(CLOCK_REALTIME, &child->mtime);
	child->mode = (S_IFDIR | 0777);

	flipBit(inodebitmap, child_inode_number);
	superblock->free_inodes_count -= 1;
	free(path_copy);
	return 0;
}



/**
 * Remove a directory.
 *
 * Implements the rmdir() system call.
 *
 * Assumptions (already verified by FUSE using getattr() calls):
 *   "path" exists and is a directory.
 *
 * Errors:
 *   ENOTEMPTY  the directory is not empty.
 *
 * @param path  path to the directory to remove.
 * @return      0 on success; -errno on error.
 */
static int a1fs_rmdir(const char *path)
{
	fs_ctx* fs = get_fs();
	a1fs_superblock* superblock = (a1fs_superblock*)(fs->image + A1FS_BLOCK_SIZE);
	a1fs_blk_t inode_table_block = superblock->inode_table;
	a1fs_inode* inodes = (a1fs_inode*)(fs->image + (A1FS_BLOCK_SIZE * inode_table_block));
	a1fs_blk_t* blockbitmap = (a1fs_blk_t*)(fs->image + (superblock->block_bitmap * A1FS_BLOCK_SIZE));
	a1fs_blk_t* inodebitmap = (a1fs_blk_t*)(fs->image + (superblock->inode_bitmap * A1FS_BLOCK_SIZE));

	//Must first check if legal (namely the dir we want to remove is empty)
	int* path_info = path_lookup(path);
	int child_inode_number = path_info[1];
	a1fs_inode* child = &inodes[child_inode_number];
	if (child->size != 0) { return -ENOTEMPTY; }

	//Now we need to get the parent inode from which we want to remove the empty dir
	int parent_inode_number = path_info[3];
	a1fs_inode* parent = &inodes[parent_inode_number];
	a1fs_extent_table* parent_extent_table = (a1fs_extent_table*)(fs->image + (parent->extent_table * A1FS_BLOCK_SIZE));

	// We want to know which index in parent to delete
	int child_index = path_info[4];

	//Scroll to first extent to change 
	int curr_extent = 0;
	int last_extent = lastExtentIndex(parent->size, parent_extent_table);

	int blocks_to_pass = (int)(child_index * sizeof(a1fs_dentry)) / A1FS_BLOCK_SIZE;
	int blocks_passed = 0;
	while (blocks_passed <= blocks_to_pass) { //Scroll to extent in which we want to remove the dir
		blocks_passed += parent_extent_table->i_extents[curr_extent].count;
		curr_extent++;
	}
	curr_extent--;
	blocks_passed -= parent_extent_table->i_extents[curr_extent].count;

	int dest_offset = (sizeof(a1fs_dentry) * child_index) - (blocks_passed * A1FS_BLOCK_SIZE);

	//Scroll to very last extent/a1fs_dentry which we can use to "fill" the hole 
	int src_offset = parent->size % A1FS_BLOCK_SIZE;
	src_offset += (parent_extent_table->i_extents[last_extent].count - 1) * A1FS_BLOCK_SIZE - sizeof(a1fs_dentry);
	if (last_extent == curr_extent && src_offset == dest_offset) { //We are removing the last dir
		parent->size -= sizeof(a1fs_dentry);
		parent->links -= 1;
		clock_gettime(CLOCK_REALTIME, &parent->mtime);
		if (parent->size % A1FS_BLOCK_SIZE == 0) { //We shrunk the extent by one block
			parent_extent_table->i_extents[last_extent].count -= 1;
			flipBit(blockbitmap, parent_extent_table->i_extents[last_extent].start + parent_extent_table->i_extents[last_extent].count);
			superblock->free_blocks_count += 1;
		}
		if (parent->size == 0) {
			flipBit(blockbitmap, parent->extent_table);
			superblock->free_blocks_count += 1;
		}
		flipBit(inodebitmap, child_inode_number);
		superblock->free_inodes_count += 1;
		return 0;
	}
	memcpy(fs->image + (parent_extent_table->i_extents[curr_extent].start * A1FS_BLOCK_SIZE) + dest_offset, fs->image + (parent_extent_table->i_extents[last_extent].start * A1FS_BLOCK_SIZE) + src_offset, sizeof(a1fs_dentry));
	parent->size -= sizeof(a1fs_dentry);
	if (parent->size % A1FS_BLOCK_SIZE == 0) { //We shrunk the extent by one block
		parent_extent_table->i_extents[last_extent].count -= 1;
		flipBit(blockbitmap, parent_extent_table->i_extents[last_extent].start + parent_extent_table->i_extents[last_extent].count);
		superblock->free_blocks_count += 1;
	}
	parent->links -= 1;
	flipBit(inodebitmap, child_inode_number);
	clock_gettime(CLOCK_REALTIME, &parent->mtime);
	return 0;
}


/**
 * Create a file.
 *
 * Implements the open()/creat() system call.
 *
 * Assumptions (already verified by FUSE using getattr() calls):
 *   "path" doesn't exist.
 *   The parent directory of "path" exists and is a directory.
 *   "path" and its components are not too long.
 *
 * Errors:
 *   ENOMEM  not enough memory (e.g. a malloc() call failed).
 *   ENOSPC  not enough free space in the file system.
 *
 * @param path  path to the file to create.
 * @param mode  file mode bits.
 * @param fi    unused.
 * @return      0 on success; -errno on error.
 */
static int a1fs_create(const char *path, mode_t mode, struct fuse_file_info *fi)
{
	(void)fi;// unused
	assert(S_ISREG(mode));
	fs_ctx* fs = get_fs();


	a1fs_superblock* superblock = (a1fs_superblock*)(fs->image + A1FS_BLOCK_SIZE);
	a1fs_blk_t inode_table_block = superblock->inode_table;
	a1fs_inode* inodes = (a1fs_inode*)(fs->image + (A1FS_BLOCK_SIZE * inode_table_block));
	a1fs_blk_t* blockbitmap = (a1fs_blk_t*)(fs->image + (superblock->block_bitmap * A1FS_BLOCK_SIZE));
	a1fs_blk_t* inodebitmap = (a1fs_blk_t*)(fs->image + (superblock->inode_bitmap * A1FS_BLOCK_SIZE));
	a1fs_blk_t blocks_needed = 0;

	//First we have to check if there are enoguh free inodes and blocks. 
	int* path_info = path_lookup(path);
	//We begin by checking if the last existing inode (to which we add more dirs) has enough space for another a1fs_dir_entry
	if ((inodes[path_info[1]].size == 0)) {
		blocks_needed++;
	}
	if ((inodes[path_info[1]].size % A1FS_BLOCK_SIZE) < sizeof(a1fs_dentry)) {
		blocks_needed++;
	}

	//Now we reconstruct the path at which path_looup exited (that is the rest of the directories we want to create)
	char* path_copy = malloc((strlen(path)) * sizeof(char));
	if (path_copy == NULL) { return -ENOMEM; }
	strcpy(path_copy, path);
	const char delimiter[2] = "/";
	char* path_left;
	path_left = strtok(path_copy, delimiter);
	for (int i = 0; i < path_info[2]; i++) { path_left = strtok(NULL, delimiter); }


	if (superblock->free_blocks_count < blocks_needed || superblock->free_inodes_count < 1) {
		free(path_copy);
		return -ENOSPC;
	}
	//------------------------------------------------------------------------------------------------------
	//Now we have the rest of a path that we need to create dirs for. 
	char* dir_to_create = strtok(path_left, delimiter);

	//Parent dir vars
	a1fs_ino_t parent_inode_number = (a1fs_ino_t)path_info[1];
	a1fs_inode* parent = &inodes[parent_inode_number];

	//Child dir vars
	int i = findZeroBit(inodebitmap, superblock->inodes_count); //is non negative thx to previous check
	printf("zero bit for inode at: %d \n", i);
	a1fs_ino_t child_inode_number = (a1fs_ino_t)(i);
	a1fs_inode* child = &inodes[child_inode_number];

	//Three cases. Case 0: parent is empty and doesnt have an extent block yet. Case 1:  we need to create a new extent, Case 2 there is enough space in current extent
	//Case 0
	if (parent->size == 0) {
		printf("case 0 \n");
		int i = findZeroBit(blockbitmap, superblock->blocks_count); //is non negative thx to previous check
		parent->extent_table = i;
		a1fs_extent_table* parent_extent_table = (a1fs_extent_table*)(fs->image + (parent->extent_table * A1FS_BLOCK_SIZE));
		flipBit(blockbitmap, i);

		//We need to set the values of the extent
		int k = findZeroBit(blockbitmap, superblock->blocks_count);
		a1fs_blk_t k2 = (a1fs_blk_t)k;
		printf("created extent table and set extent start to %d \n", k2);
		parent_extent_table->i_extents[0].start = k2;
		parent_extent_table->i_extents[0].count = 1;
		flipBit(blockbitmap, k2);
		superblock->free_blocks_count -= 2;

	}

	int space_left = (int)(A1FS_BLOCK_SIZE - (parent->size % A1FS_BLOCK_SIZE));
	int space_used = (int)A1FS_BLOCK_SIZE - space_left;
	int extent_offset = 0;
	if (space_left == A1FS_BLOCK_SIZE && parent->size != 0) { //Case 1
		printf("-----------case 1-----------\n");
		a1fs_extent_table* parent_extent_table = (a1fs_extent_table*)(fs->image + (parent->extent_table * A1FS_BLOCK_SIZE));
		int last_extent = lastExtentIndex((int)(parent->size), parent_extent_table);
		int last_block = parent_extent_table->i_extents[last_extent].start + parent_extent_table->i_extents[last_extent].count - 1;
		//We try to extend last last extent by one to keep fragmentation low
		int k = checkZeroBit(blockbitmap, last_block + 1, superblock->blocks_count);
		if (k == 0) {
			parent_extent_table->i_extents[last_extent].count += 1;
			flipBit(blockbitmap, last_block + 1);
		}
		else { // Need to create new extent
			if (last_extent + 1 >= 512) {
				free(path_copy);
				return -ENOSPC;
			}
			int i = findZeroBit(blockbitmap, superblock->blocks_count);
			printf("find zero bit in case 1: %d \n", i);
			parent_extent_table->i_extents[last_extent + 1].start = i;
			printf("extent start in case 1: %d \n", parent_extent_table->i_extents[last_extent + 1].start);
			parent_extent_table->i_extents[last_extent + 1].count = 1;
			flipBit(blockbitmap, i);
			extent_offset += 1;
		}

	}
	//Go to last position in extent table (using size of inode)
	a1fs_extent_table* parent_extent_table = (a1fs_extent_table*)(fs->image + (parent->extent_table * A1FS_BLOCK_SIZE));
	printf("case 2 \n");
	int last_extent = lastExtentIndex((int)(parent->size), parent_extent_table) + extent_offset;
	printf("last extent: %d \n", last_extent);
	printf("extent offset: %d \n", extent_offset);
	int offset = (int)(space_used);
	int last_block = (int)(parent_extent_table->i_extents[last_extent].start + parent_extent_table->i_extents[last_extent].count) - 1;
	printf("last extent start: %d \n", parent_extent_table->i_extents[last_extent].start);
	printf("offset: %d \n", offset);
	printf("last block: %d \n", last_block);

	//Now we can add the dir entry 
	a1fs_dentry* dentry = fs->image + (last_block * A1FS_BLOCK_SIZE) + offset;
	strncpy(dentry->name, dir_to_create, A1FS_NAME_MAX);
	dentry->ino = child_inode_number;

	//Now we update the parent's info
	parent->size = parent->size + sizeof(a1fs_dentry);
	parent->links = parent->links + 1;
	clock_gettime(CLOCK_REALTIME, &parent->mtime);

	//Now we update the child's info
	child->size = 0;
	child->links = 1;
	clock_gettime(CLOCK_REALTIME, &child->mtime);
	child->mode = mode;

	flipBit(inodebitmap, child_inode_number);
	superblock->free_inodes_count -= 1;
	free(path_copy);
	return 0;
}

/**
 * Remove a file.
 *
 * Implements the unlink() system call.
 *
 * Assumptions (already verified by FUSE using getattr() calls):
 *   "path" exists and is a file.
 *
 * Errors: none
 *
 * @param path  path to the file to remove.
 * @return      0 on success; -errno on error.
 */
static int a1fs_unlink(const char *path)
{
	fs_ctx* fs = get_fs();
	a1fs_superblock* superblock = (a1fs_superblock*)(fs->image + A1FS_BLOCK_SIZE);
	a1fs_blk_t inode_table_block = superblock->inode_table;
	a1fs_inode* inodes = (a1fs_inode*)(fs->image + (A1FS_BLOCK_SIZE * inode_table_block));
	a1fs_blk_t* blockbitmap = (a1fs_blk_t*)(fs->image + (superblock->block_bitmap * A1FS_BLOCK_SIZE));
	a1fs_blk_t* inodebitmap = (a1fs_blk_t*)(fs->image + (superblock->inode_bitmap * A1FS_BLOCK_SIZE));

	//Must first check if legal (namely the dir we want to remove is empty)
	int* path_info = path_lookup(path);
	int child_inode_number = path_info[1];
	a1fs_inode* child = &inodes[child_inode_number];

	//Now we need to get the parent inode from which we want to remove the empty dir
	int parent_inode_number = path_info[3];
	a1fs_inode* parent = &inodes[parent_inode_number];
	a1fs_extent_table* parent_extent_table = (a1fs_extent_table*)(fs->image + (parent->extent_table * A1FS_BLOCK_SIZE));

	child->links -= 1;
	if (child->links != 0) {
		return 0;
	}
	//NEED TO TRUNCATE *FREE table + extents
	//******************************************************************** NEED SOMEHTING HERE, NEED SOMEHTING HERE,  NEED SOMEHTING HERE,  NEED SOMEHTING HERE,  NEED SOMEHTING HERE,  NEED SOMEHTING HERE,  NEED SOMEHTING HERE,  NEED SOMEHTING HERE, 
    free_extends(child, 0);

	// We want to know which index in parent to delete
	int child_index = path_info[4];

	//Scroll to first extent to change 
	int curr_extent = 0;
	int last_extent = lastExtentIndex(parent->size, parent_extent_table);

	int blocks_to_pass = (int)(child_index * sizeof(a1fs_dentry)) / A1FS_BLOCK_SIZE;
	int blocks_passed = 0;
	while (blocks_passed <= blocks_to_pass) { //Scroll to extent in which we want to remove the dir
		blocks_passed += parent_extent_table->i_extents[curr_extent].count;
		curr_extent++;
	}
	curr_extent--;
	blocks_passed -= parent_extent_table->i_extents[curr_extent].count;

	int dest_offset = (sizeof(a1fs_dentry) * child_index) - (blocks_passed * A1FS_BLOCK_SIZE);

	//Scroll to very last extent/a1fs_dentry which we can use to "fill" the hole 
	int src_offset = parent->size % A1FS_BLOCK_SIZE;
	src_offset += (parent_extent_table->i_extents[last_extent].count - 1) * A1FS_BLOCK_SIZE - sizeof(a1fs_dentry);
	if (last_extent == curr_extent && src_offset == dest_offset) { //We are removing the last dir
		parent->size -= sizeof(a1fs_dentry);
		parent->links -= 1;
		clock_gettime(CLOCK_REALTIME, &parent->mtime);
		if (parent->size % A1FS_BLOCK_SIZE == 0) { //We shrunk the extent by one block
			parent_extent_table->i_extents[last_extent].count -= 1;
			flipBit(blockbitmap, parent_extent_table->i_extents[last_extent].start + parent_extent_table->i_extents[last_extent].count);
			superblock->free_blocks_count += 1;
		}
		if (parent->size == 0) {
			flipBit(blockbitmap, parent->extent_table);
			superblock->free_blocks_count += 1;
		}
		flipBit(inodebitmap, child_inode_number);
		superblock->free_inodes_count += 1;
		return 0;
	}
	memcpy(fs->image + (parent_extent_table->i_extents[curr_extent].start * A1FS_BLOCK_SIZE) + dest_offset, fs->image + (parent_extent_table->i_extents[last_extent].start * A1FS_BLOCK_SIZE) + src_offset, sizeof(a1fs_dentry));
	parent->size -= sizeof(a1fs_dentry);
	if (parent->size % A1FS_BLOCK_SIZE == 0) { //We shrunk the extent by one block
		parent_extent_table->i_extents[last_extent].count -= 1;
		flipBit(blockbitmap, parent_extent_table->i_extents[last_extent].start + parent_extent_table->i_extents[last_extent].count);
		superblock->free_blocks_count += 1;
	}
	flipBit(inodebitmap, child_inode_number);
	clock_gettime(CLOCK_REALTIME, &parent->mtime);
	return 0;
}


/**
 * Change the modification time of a file or directory.
 *
 * Implements the utimensat() system call. See "man 2 utimensat" for details.
 *
 * NOTE: You only need to implement the setting of modification time (mtime).
 *       Timestamp modifications are not recursive. 
 *
 * Assumptions (already verified by FUSE using getattr() calls):
 *   "path" exists.
 *
 * Errors: none
 *
 * @param path   path to the file or directory.
 * @param times  timestamps array. See "man 2 utimensat" for details.
 * @return       0 on success; -errno on failure.
 */
static int a1fs_utimens(const char *path, const struct timespec times[2])
{
	fs_ctx *fs = get_fs();

	//TODO: update the modification timestamp (mtime) in the inode for given
	// path with either the time passed as argument or the current time,
	// according to the utimensat man page
	int* path_info = path_lookup(path);
	a1fs_superblock* superblock = (a1fs_superblock*)(fs->image + A1FS_BLOCK_SIZE);
	a1fs_blk_t inode_table_block = superblock->inode_table;
	a1fs_inode* inodes = (a1fs_inode*)(fs->image + (A1FS_BLOCK_SIZE * inode_table_block));
	int child_inode_number = path_info[1];
	a1fs_inode* child = &inodes[child_inode_number];
	if (times == NULL) {
		clock_gettime(CLOCK_REALTIME, &child->mtime);
	}
	child->mtime = times[1];

	return 0;
}

/**
 * Change the size of a file.
 *
 * Implements the truncate() system call. Supports both extending and shrinking.
 * If the file is extended, the new uninitialized range at the end must be
 * filled with zeros.
 *
 * Assumptions (already verified by FUSE using getattr() calls):
 *   "path" exists and is a file.
 *
 * Errors:
 *   ENOMEM  not enough memory (e.g. a malloc() call failed).
 *   ENOSPC  not enough free space in the file system.
 *
 * @param path  path to the file to set the size.
 * @param size  new file size in bytes.
 * @return      0 on success; -errno on error.
 */
static int a1fs_truncate(const char *path, off_t size)
{
	//Getting vars
	fs_ctx *fs = get_fs();
	a1fs_superblock* superblock = (a1fs_superblock*)(fs->image + A1FS_BLOCK_SIZE);
	a1fs_blk_t inode_table_block = superblock->inode_table;
	a1fs_inode* inodes = (a1fs_inode*)(fs->image + (A1FS_BLOCK_SIZE * inode_table_block));
	a1fs_blk_t* blockbitmap = (a1fs_blk_t*)(fs->image + (superblock->block_bitmap * A1FS_BLOCK_SIZE));

	//Get file we want to truncate
	int* path_info = path_lookup(path);
	int child_inode_number = path_info[1];
	a1fs_inode* file = &inodes[child_inode_number];
	a1fs_extent_table* file_extent_table = (a1fs_extent_table*)(fs->image + (file->extent_table * A1FS_BLOCK_SIZE));

	//Two cases: 1. we are making the file smaller, 2 we are making the file bigger
	//Case 1 
	if (size < (int)file->size) {
		free_extends(file, size);
	}
	//Case 2
	if (size > (int)file->size) {
        int blocks_needed = 0;


        if ((inodes[path_info[1]].size == 0)) {
            printf("case 0 \n");
            if(superblock->free_blocks_count == 0) {return -ENOSPC;}
            superblock->free_blocks_count -= 1;

            int i = findZeroBit(blockbitmap, superblock->blocks_count); //is non negative thx to previous check
            file->extent_table = i;
            file_extent_table = (a1fs_extent_table*)(fs->image + (file->extent_table * A1FS_BLOCK_SIZE));
            flipBit(blockbitmap, i);

        }

        if ((inodes[path_info[1]].size == 0)) { // if the file size is 0
            int remaining_size = size;
            int extent_index = 0;
            while(remaining_size > 0){ //creates new extents
                int* extent_info = longest_extent(blockbitmap, superblock->blocks_count); // get the longest available extent
                file_extent_table->i_extents[extent_index].start = extent_info[0];
                int extent_count;
                if((int) extent_info[1]*A1FS_BLOCK_SIZE >= (int) size){ // if the extent is big enough to take the whole data
                    extent_count = (size)/A1FS_BLOCK_SIZE + 1;
                }
                else{
                    extent_count = extent_info[1];
                }
                for(int x = 0; x < extent_count; x++){flipBit(blockbitmap, x + extent_info[0]);}
                file_extent_table->i_extents[extent_index].count = (int) (extent_count);
                remaining_size -= extent_count*A1FS_BLOCK_SIZE;
                extent_index++;
                blocks_needed += extent_count;
            }
            inodes[path_info[1]].size = size;
        }
        else{ // if file is not empty
            int remaining_size = size - (int)file->size;

            int last_ext = lastExtentIndex(inodes[path_info[1]].size, file_extent_table);

            int if_last_one_is_full = 0; // a flag to show if the last extend cannot be extended
            while(remaining_size > 0){
                int last_block = file_extent_table->i_extents[last_ext].start + file_extent_table->i_extents[last_ext].count;
                if ((checkZeroBit(blockbitmap, last_block + 1, superblock->blocks_count) == 0) &&
                    if_last_one_is_full == 0) {
                    file_extent_table->i_extents[last_ext].count += 1;
                    remaining_size -= A1FS_BLOCK_SIZE;
                    blocks_needed += 1;
                    flipBit(blockbitmap, last_block + 1);
                }
                else{
                    if_last_one_is_full = 1;
                    last_ext++;
                    int *extent_info = longest_extent(blockbitmap,
                                                      superblock->blocks_count); // get the longest available extent
                    file_extent_table->i_extents[last_ext].start = extent_info[0];
                    int extent_count;
                    if (extent_info[1] * A1FS_BLOCK_SIZE >=
                        remaining_size) { // if the extent is big enough to take the whole data
                        extent_count = ((int) remaining_size) / ((int) A1FS_BLOCK_SIZE) + 1;
                    } else {
                        extent_count = extent_info[1];
                    }
                    for(int x = 0; x < extent_count; x++){flipBit(blockbitmap, x + extent_info[0]);}
                    file_extent_table->i_extents[last_ext].count = (int) (extent_count);
                    remaining_size -= extent_count * A1FS_BLOCK_SIZE;
                    blocks_needed += extent_count;
                }
                printf("case 5 \n");
            }
            inodes[path_info[1]].size = size;
                    }



	}

	return 0;
}


/**
 * Read data from a file.
 *
 * Implements the pread() system call. Must return exactly the number of bytes
 * requested except on EOF (end of file). Reads from file ranges that have not
 * been written to must return ranges filled with zeros. You can assume that the
 * byte range from offset to offset + size is contained within a single block.
 *
 * Assumptions (already verified by FUSE using getattr() calls):
 *   "path" exists and is a file.
 *
 * Errors: none
 *
 * @param path    path to the file to read from.
 * @param buf     pointer to the buffer that receives the data.
 * @param size    buffer size (number of bytes requested).
 * @param offset  offset from the beginning of the file to read from.
 * @param fi      unused.
 * @return        number of bytes read on success; 0 if offset is beyond EOF;
 *                -errno on error.
 */
static int a1fs_read(const char *path, char *buf, size_t size, off_t offset,
                     struct fuse_file_info *fi)
{
	(void)fi;// unused
	fs_ctx *fs = get_fs();

	//TODO: read data from the file at given offset into the buffer
    a1fs_superblock* superblock = (a1fs_superblock*)(fs->image + A1FS_BLOCK_SIZE);
    a1fs_blk_t inode_table_block = superblock->inode_table;
    a1fs_inode* inodes = (a1fs_inode*)(fs->image + (A1FS_BLOCK_SIZE * inode_table_block));

    //Get file we want to truncate
    int* path_info = path_lookup(path);
    int child_inode_number = path_info[1];
    a1fs_inode* file = &inodes[child_inode_number];
    a1fs_extent_table* file_extent_table = (a1fs_extent_table*)(fs->image + (file->extent_table * A1FS_BLOCK_SIZE));
    int last_ext = lastExtentIndex(file->size, file_extent_table);
    if(offset >= (int) file->size){return 0;}
    int curr_size = 0;
    int start_found = 0;
    int amount_copied = 0;
    for(int i = 0; i <= last_ext; i++){
        if(start_found == 1){
            memcpy(&buf[amount_copied],file_extent_table->i_extents[i].start * A1FS_BLOCK_SIZE + fs->image, file_extent_table->i_extents[i].count * A1FS_BLOCK_SIZE);
            amount_copied += (file_extent_table->i_extents[i].count * A1FS_BLOCK_SIZE);
        }
        if(offset >= curr_size && offset < curr_size+file_extent_table->i_extents[i].count * A1FS_BLOCK_SIZE){
            memcpy(buf,file_extent_table->i_extents[i].start * A1FS_BLOCK_SIZE + offset + fs->image, file_extent_table->i_extents[i].count * A1FS_BLOCK_SIZE);
            amount_copied += ((file_extent_table->i_extents[i].count * A1FS_BLOCK_SIZE)-(offset));
            start_found = 1;
        }
        curr_size += file_extent_table->i_extents[i].count * A1FS_BLOCK_SIZE;
    }
	return amount_copied;
}

/**
 * Write data to a file.
 *
 * Implements the pwrite() system call. Must return exactly the number of bytes
 * requested except on error. If the offset is beyond EOF (end of file), the
 * file must be extended. If the write creates a "hole" of uninitialized data,
 * the new uninitialized range must filled with zeros. You can assume that the
 * byte range from offset to offset + size is contained within a single block.
 *
 * Assumptions (already verified by FUSE using getattr() calls):
 *   "path" exists and is a file.
 *
 * Errors:
 *   ENOMEM  not enough memory (e.g. a malloc() call failed).
 *   ENOSPC  not enough free space in the file system.
 *   ENOSPC  too many extents (a1fs only needs to support 512 extents per file)
 *
 * @param path    path to the file to write to.
 * @param buf     pointer to the buffer containing the data.
 * @param size    buffer size (number of bytes requested).
 * @param offset  offset from the beginning of the file to write to.
 * @param fi      unused.
 * @return        number of bytes written on success; -errno on error.
 */
static int a1fs_write(const char *path, const char *buf, size_t size,
                      off_t offset, struct fuse_file_info *fi)
{
	(void)fi;// unused
	fs_ctx *fs = get_fs();

	//TODO: write data from the buffer into the file at given offset, possibly
	// "zeroing out" the uninitialized range
    int* path_info = path_lookup(path);
    //a1fs_ino_t inode_num = (a1fs_ino_t) path_info[1];
    a1fs_superblock* superblock = (a1fs_superblock*)(fs->image + A1FS_BLOCK_SIZE);
    a1fs_blk_t inode_table_block = superblock->inode_table;
    a1fs_inode* inodes = (a1fs_inode*)(fs->image + (A1FS_BLOCK_SIZE * inode_table_block));
    a1fs_blk_t* blockbitmap = (a1fs_blk_t*)(fs->image + (superblock->block_bitmap * A1FS_BLOCK_SIZE));
    //a1fs_blk_t* inodebitmap = (a1fs_blk_t*)(image + (superblock->inode_bitmap * A1FS_BLOCK_SIZE));
    a1fs_blk_t* data_blocks = (a1fs_blk_t*)(fs->image + (superblock->data_blocks * A1FS_BLOCK_SIZE));
    a1fs_blk_t blocks_needed = 0;

    int num_bytes_written = 0;
    int inode_number = path_info[1];
    a1fs_inode* this_ino = &inodes[inode_number];
    a1fs_extent_table* extent_table = (a1fs_extent_table*)(fs->image + (this_ino->extent_table * A1FS_BLOCK_SIZE));

    if ((inodes[path_info[1]].size == 0)) {
        printf("case 0 \n");
        if(superblock->free_blocks_count == 0) {return -ENOSPC;}
        superblock->free_blocks_count -= 1;

        int i = findZeroBit(blockbitmap, superblock->blocks_count); //is non negative thx to previous check
        this_ino->extent_table = i;
        extent_table = (a1fs_extent_table*)(fs->image + (this_ino->extent_table * A1FS_BLOCK_SIZE));
        flipBit(blockbitmap, i);

    }

    if ((inodes[path_info[1]].size == 0)) { // if the file size is 0
        char *words = create_string(buf, offset, 0);
        num_bytes_written = strlen(words);
        int remaining_size = offset+size;
        int extent_index = 0;
        int copy_start = 0;
        while(remaining_size > 0){ //creates new extents
            int* extent_info = longest_extent(blockbitmap, superblock->blocks_count); // get the longest available extent
            extent_table->i_extents[extent_index].start = extent_info[0];
            int extent_count;
            if((int) extent_info[1]*A1FS_BLOCK_SIZE >= (int) offset+size){ // if the extent is big enough to take the whole data
                extent_count = (offset+size)/A1FS_BLOCK_SIZE + 1;
            }
            else{
                extent_count = extent_info[1];
            }
            for(int x = 0; x < extent_count; x++){flipBit(blockbitmap, x + extent_info[0]);}
            extent_table->i_extents[extent_index].count = (int) (extent_count);
            memcpy(fs->image + (extent_info[0] * A1FS_BLOCK_SIZE), &words[copy_start], extent_count * A1FS_BLOCK_SIZE);
            copy_start += extent_count * A1FS_BLOCK_SIZE;
            remaining_size -= extent_count*A1FS_BLOCK_SIZE;
            extent_index++;
            blocks_needed += extent_count;
        }
        printf("case 1 \n");
        inodes[path_info[1]].size = offset+size;
    }
    else{ // if file is not empty
        printf("case 3 \n");
        char *words = create_string(buf, offset, (int) inodes[path_info[1]].size); //creates a string filled with zeros
        num_bytes_written = strlen(words);
        // in the beginning if there is a gap between previous end of the file and the offset
        int new_offset;
        int remaining_size;
        if(offset > (int) inodes[path_info[1]].size){
            new_offset = (int) inodes[path_info[1]].size;
            remaining_size = (offset - new_offset) + size;
        }
        else{
            new_offset = offset;
            remaining_size = size;
        } // set new offset equal to the lesser of the current end of the file and the offset
        // so that the new starting point will be in range of the extents


        int last_ext = lastExtentIndex(inodes[path_info[1]].size, extent_table);
        int ext_index = 0;
        int curr_extent_size = 0;

        while((curr_extent_size + (extent_table->i_extents[ext_index].count) * A1FS_BLOCK_SIZE ) < new_offset){ // find the extent to start writing
            ext_index++;
            curr_extent_size += (extent_table->i_extents[ext_index].count) * A1FS_BLOCK_SIZE;
        }

        int if_last_one_is_full = 0; // a flag to show if the last extend cannot be extended
        int copy_index = 0; // index of the words to start copying
        while(remaining_size > 0){
            if(ext_index <= last_ext) {
                memcpy(fs->image + (extent_table->i_extents[ext_index].start * A1FS_BLOCK_SIZE) + (new_offset - curr_extent_size), &words[copy_index], extent_table->i_extents[ext_index].count * A1FS_BLOCK_SIZE);
                if(extent_table->i_extents[ext_index].count * A1FS_BLOCK_SIZE < new_offset - curr_extent_size + remaining_size) {
                    ext_index++;
                }
                int amount_copied = extent_table->i_extents[ext_index].count * A1FS_BLOCK_SIZE - (new_offset - curr_extent_size);
                copy_index += (amount_copied);
                remaining_size -= amount_copied;
                printf("case 4 \n");
            }
            else{
                int last_block = extent_table->i_extents[last_ext].start + extent_table->i_extents[last_ext].count;
                if ((checkZeroBit(blockbitmap, last_block + 1, superblock->blocks_count) == 0) &&
                    if_last_one_is_full == 0) {
                    extent_table->i_extents[last_ext].count += 1;
                    remaining_size -= A1FS_BLOCK_SIZE;
                    memcpy(fs->image + ((extent_table->i_extents[ext_index].start+extent_table->i_extents[last_ext].count-1) * A1FS_BLOCK_SIZE), &words[copy_index],  A1FS_BLOCK_SIZE);
                    copy_index += A1FS_BLOCK_SIZE;
                    blocks_needed += 1;
                    flipBit(blockbitmap, last_block + 1);
                }
                else{
                    if_last_one_is_full = 1;
                    last_ext++;
                    int *extent_info = longest_extent(blockbitmap,
                                                      superblock->blocks_count); // get the longest available extent
                    extent_table->i_extents[last_ext].start = extent_info[0];
                    int extent_count;
                    if (extent_info[1] * A1FS_BLOCK_SIZE >=
                        (offset + size)) { // if the extent is big enough to take the whole data
                        extent_count = ((int) (offset + size)) / ((int) A1FS_BLOCK_SIZE) + 1;
                    } else {
                        extent_count = extent_info[1];
                    }
                    for(int x = 0; x < extent_count; x++){flipBit(blockbitmap, x + extent_info[0]);}
                    extent_table->i_extents[last_ext].count = (int) (extent_count);
                    remaining_size -= extent_count * A1FS_BLOCK_SIZE;
                    memcpy(fs->image + (extent_info[0] * A1FS_BLOCK_SIZE), &words[copy_index], extent_count * A1FS_BLOCK_SIZE);
                    copy_index += extent_count * A1FS_BLOCK_SIZE;
                    blocks_needed += extent_count;
                }
                printf("case 5 \n");
            }
            if(inodes[path_info[1]].size < offset+size){
                inodes[path_info[1]].size = offset+size;
            }
        }
    }
    superblock->free_blocks_count -= blocks_needed;
    return num_bytes_written;
}


static struct fuse_operations a1fs_ops = {
	.destroy  = a1fs_destroy,
	.statfs   = a1fs_statfs,
	.getattr  = a1fs_getattr,
	.readdir  = a1fs_readdir,
	.mkdir    = a1fs_mkdir,
	.rmdir    = a1fs_rmdir,
	.create   = a1fs_create,
	.unlink   = a1fs_unlink,
	.utimens  = a1fs_utimens,
	.truncate = a1fs_truncate,
	.read     = a1fs_read,
	.write    = a1fs_write,
};

int main(int argc, char *argv[])
{
	a1fs_opts opts = {0};// defaults are all 0
	struct fuse_args args = FUSE_ARGS_INIT(argc, argv);
	if (!a1fs_opt_parse(&args, &opts)) return 1;

	fs_ctx fs = {0};
	if (!a1fs_init(&fs, &opts)) {
		fprintf(stderr, "Failed to mount the file system\n");
		return 1;
	}

	return fuse_main(args.argc, args.argv, &a1fs_ops, &fs);
}






