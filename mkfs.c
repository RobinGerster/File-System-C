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
 * CSC369 Assignment 1 - a1fs formatting tool.
 */

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <unistd.h>

#include "a1fs.h"
#include "map.h"
#include "time.h"


/** Command line options. */
typedef struct mkfs_opts {
	/** File system image file path. */
	const char *img_path;
	/** Number of inodes. */
	size_t n_inodes;

	/** Print help and exit. */
	bool help;
	/** Overwrite existing file system. */
	bool force;
	/** Zero out image contents. */
	bool zero;

} mkfs_opts;

static const char *help_str = "\
Usage: %s options image\n\
\n\
Format the image file into a1fs file system. The file must exist and\n\
its size must be a multiple of a1fs block size - %zu bytes.\n\
\n\
Options:\n\
    -i num  number of inodes; required argument\n\
    -h      print help and exit\n\
    -f      force format - overwrite existing a1fs file system\n\
    -z      zero out image contents\n\
";

static void print_help(FILE *f, const char *progname)
{
	fprintf(f, help_str, progname, A1FS_BLOCK_SIZE);
}


static bool parse_args(int argc, char *argv[], mkfs_opts *opts)
{
	char o;
	while ((o = getopt(argc, argv, "i:hfvz")) != -1) {
		switch (o) {
			case 'i': opts->n_inodes = strtoul(optarg, NULL, 10); break;

			case 'h': opts->help  = true; return true;// skip other arguments
			case 'f': opts->force = true; break;
			case 'z': opts->zero  = true; break;

			case '?': return false;
			default : assert(false);
		}
	}

	if (optind >= argc) {
		fprintf(stderr, "Missing image path\n");
		return false;
	}
	opts->img_path = argv[optind];

	if (opts->n_inodes == 0) {
		fprintf(stderr, "Missing or invalid number of inodes\n");
		return false;
	}
	return true;
}


/** Determine if the image has already been formatted into a1fs. */
static bool a1fs_is_present(void *image)
{
	//TODO: check if the image already contains a valid a1fs superblock
	a1fs_superblock* superblock = (a1fs_superblock*)(image + A1FS_BLOCK_SIZE);
	if (superblock->magic == A1FS_MAGIC) {
		return true;
	}
	else {
		return false;
	}
}



//A function to flip a bit. flipt the kth bit in bitmap.
void  flipBit(a1fs_blk_t* bitmap,  int k)
{
	int i = k / 32;        //gives the corresponding index in the array A
	int pos = k % 32;      //gives the corresponding bit position in A[i]

	unsigned int flag = 1;   // flag = 0000.....00001
	flag = flag << pos;      // flag = 0000...010...000   (shifted k positions)

	bitmap[i] = bitmap[i] ^ flag;      // Flip the bit at the k-th position in A[i]
}


//find the first 0 in the bitmap where totalBlocks = superblock->block_count. Returns the block number of the corresponding free block.
int findZeroBit(a1fs_blk_t* bitmap, a1fs_blk_t totalBlocks) {
	for (a1fs_blk_t i = 0; i < totalBlocks/32+1; i++) {
		if (bitmap[i] != 2147483647) {
			for (a1fs_blk_t k = 0; k < 32; k++) {//gives the corresponding bit position in A[i]
				unsigned int flag = 1;   // flag = 0000.....00001
				flag = flag << k;      // flag = 0000...010...000   (shifted k positions)
				if ((bitmap[i] & flag) == 0) { //Then the bit at pos k is free
					a1fs_blk_t free = (i * 32 + k);
					if (free >= totalBlocks) {
						return -1;
					}
					return (int)(i * 32 + k);
				}
			}
		}
	}
	return -1;
}




//DEBUG STUFF -----------------------------------------------------------------------------------------------------------------------
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
    for (a1fs_blk_t i = 0; i < totalBlocks/32+1; i++) {
        if (bitmap[i] != 2147483647) {
            for (a1fs_blk_t k = 0; k < 32; k++) {//gives the corresponding bit position in A[i]
                unsigned int flag = 1;   // flag = 0000.....00001
                flag = flag << k;      // flag = 0000...010...000   (shifted k positions)
                if ((bitmap[i] & flag) == 0 && i * 32 + k < totalBlocks) { //Then the bit at pos k is free
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
    if(longest == current_length){ // if this was the longest extent
        best_start = current_best_start; // set best start to the current best start
    }
    if(longest > 0){
        toReturn[0] = best_start;
        toReturn[1] = longest;
        return toReturn;
    }
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
        memset(new_buf,'0',zeros);
        memcpy(&new_buf[zeros], buf, total_length);
        return new_buf;
    }
}

int* path_lookup(const char* path, void* image) {
	int* toReturn = malloc(sizeof(int) * 3);
	toReturn[0] = 0;
	toReturn[1] = 0;
	toReturn[2] = 0;
	toReturn[3] = 0;
	toReturn[4] = 0;
	// Get the array of inodes (inode table)
	a1fs_superblock* superblock = (a1fs_superblock*)(image + A1FS_BLOCK_SIZE);
	a1fs_blk_t inode_table_block = superblock->inode_table;
	a1fs_inode* inodes = (a1fs_inode*)(image + (A1FS_BLOCK_SIZE * inode_table_block));
	int inode_number = 0;
	int parent_number = 0;

	printf("\n");
	printf("\n");

	//Get first token in path
	char* path_copy = malloc((strlen(path)) * sizeof(char));
	strcpy(path_copy, path);
	char* dir_entry_name;
	const char delimiter[2] = "/";
	dir_entry_name = strtok(path_copy, delimiter); // Get first element in path

	//Go to last token in path
	while (dir_entry_name != NULL) {
		//printf("looking for: %s \n", dir_entry_name);
		//Check if inode is a valid dir (if not switch index 0 of toReturn to -1 for ENOTDIR)
		if (((&inodes[inode_number])->mode & S_IFDIR) != S_IFDIR) {
			printf("parent identified as regular file \n");
			toReturn[0] = -1;
			toReturn[1] = inode_number;
			toReturn[3] = parent_number;
			free(path_copy);
			return toReturn;
		}
		//Some vars we need to keep track of
		int dir_entry_count = inodes[inode_number].size / sizeof(a1fs_dentry); // How many dir entries does the current inode have
		a1fs_extent_table* extent_table = (a1fs_extent_table*)(image + (inodes[inode_number].extent_table * A1FS_BLOCK_SIZE)); //The array of extents corresponding to the inode
		int extent_number = 0;
		a1fs_dentry* dentry = (a1fs_dentry*)(image + (A1FS_BLOCK_SIZE * extent_table->i_extents[0].start));
		//Loop over all dir entries
		int error = -2;  //Want to return this for ENOENT
		int k = 0; //The index into a dir_entry block (this reflects changes when going to a new extent)
		for (int i = 0; i < dir_entry_count; i++) {
			//printf("comparing %s ", dentry[k].name);
			if (strcmp(dentry[k].name, dir_entry_name) == 0) { //Is the dentry.name equal to the current token?
				parent_number = inode_number;
				toReturn[4] = i;
				inode_number = (int)dentry[k].ino;
				extent_number = 0;
				error = 0;
				printf("found to be equal \n");
				break;
			}
			k++;
			if (k * sizeof(a1fs_dentry) == extent_table->i_extents[extent_number].count * A1FS_BLOCK_SIZE) {  //Is the next dir entry in a new extent?
				extent_number++;
				dentry = (a1fs_dentry*)(image + (A1FS_BLOCK_SIZE * extent_table->i_extents[extent_number].start));
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
	printf("\n");
	return toReturn;
}

int mkdir2(const char* path, void* image)
{
	printf("\n");
	printf("mkdir begun \n");
	a1fs_superblock* superblock = (a1fs_superblock*)(image + A1FS_BLOCK_SIZE);
	a1fs_blk_t inode_table_block = superblock->inode_table;
	a1fs_inode* inodes = (a1fs_inode*)(image + (A1FS_BLOCK_SIZE * inode_table_block));
	a1fs_blk_t* blockbitmap = (a1fs_blk_t*)(image + (superblock->block_bitmap * A1FS_BLOCK_SIZE));
	a1fs_blk_t* inodebitmap = (a1fs_blk_t*)(image + (superblock->inode_bitmap * A1FS_BLOCK_SIZE));
	a1fs_blk_t blocks_needed = 0;

	//First we have to check if there are enoguh free inodes and blocks.
	int* path_info = path_lookup(path, image);
	//We begin by checking if the last existing inode (to which we add more dirs) has enough space for another a1fs_dir_entry
	if ((inodes[path_info[1]].size == 0)) {
		blocks_needed++;
	}
	if ((inodes[path_info[1]].size % A1FS_BLOCK_SIZE) < sizeof(a1fs_dentry)) {
		blocks_needed++;
	}

	//Now we reconstruct the path at which path_looup exited (that is the rest of the directories we want to create)
	char* path_copy = malloc((strlen(path)) * sizeof(char));
	if (path_copy == NULL) { return -1; }
	strcpy(path_copy, path);
	const char delimiter[2] = "/";
	char* path_left;
	path_left = strtok(path_copy, delimiter);
	for (int i = 0; i < path_info[2]; i++) { path_left = strtok(NULL, delimiter); }


	if (superblock->free_blocks_count < blocks_needed || superblock->free_inodes_count < 1) {
		free(path_copy);
		return -1;
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
		a1fs_extent_table* parent_extent_table = (a1fs_extent_table*)(image + (parent->extent_table * A1FS_BLOCK_SIZE));
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
		a1fs_extent_table* parent_extent_table = (a1fs_extent_table*)(image + (parent->extent_table * A1FS_BLOCK_SIZE));
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
				return -1;
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
	a1fs_extent_table* parent_extent_table = (a1fs_extent_table*)(image + (parent->extent_table * A1FS_BLOCK_SIZE));
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
	a1fs_dentry* dentry = image + (last_block * A1FS_BLOCK_SIZE) + offset;
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


int create2(const char *path, mode_t mode, void *img)
{
    assert(S_ISREG(mode));
    //TODO: create a file at given path with given mode
    char* path_copy = malloc((strlen(path)) * sizeof(char));
    if (path_copy == NULL) { return -1; } // change it to -ENOMEM
    strcpy(path_copy, path);
    const char delimiter[2] = "/";
    char* file_name;
    file_name = strtok(path_copy, delimiter);
    char file_name_actual[A1FS_NAME_MAX];
    while(file_name != NULL){
        strncpy(file_name_actual, file_name, A1FS_NAME_MAX);
        file_name = strtok(NULL, delimiter);
    }
    printf("%s \n", file_name_actual);

    int* path_info = path_lookup(path, img);
    a1fs_ino_t inode_num = (a1fs_ino_t) path_info[1];
    a1fs_superblock* superblock = (a1fs_superblock*)(img + A1FS_BLOCK_SIZE);
    a1fs_blk_t inode_table_block = superblock->inode_table;
    a1fs_inode* inodes = (a1fs_inode*)(img + (A1FS_BLOCK_SIZE * inode_table_block));
    a1fs_inode* parent = &inodes[inode_num];
    a1fs_blk_t* blockbitmap = (a1fs_blk_t*)(img + (superblock->block_bitmap * A1FS_BLOCK_SIZE));
    a1fs_blk_t* inodebitmap = (a1fs_blk_t*)(img + (superblock->inode_bitmap * A1FS_BLOCK_SIZE));

    //Child dir vars
    int i = findZeroBit(inodebitmap, superblock->inodes_count); //is non negative thx to previous check
    a1fs_ino_t child_inode_number = (a1fs_ino_t)(i);
    a1fs_inode* child = &inodes[child_inode_number];

    if (parent->size == 0) {
        printf("case 0 \n");
        int i = findZeroBit(blockbitmap, superblock->blocks_count); //is non negative thx to previous check
        parent->extent_table = i;
        a1fs_extent_table* parent_extent_table = (a1fs_extent_table*)(img + (parent->extent_table * A1FS_BLOCK_SIZE));
        flipBit(blockbitmap, i); //update the bitmap

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
    if (space_left >= (int)sizeof(a1fs_dentry)) { //Case 1
        //Go to last position in extent table (using size of inode)
        printf("case 2 \n");
        a1fs_extent_table *parent_extent_table = (a1fs_extent_table *) (img +
                                                                        (parent->extent_table * A1FS_BLOCK_SIZE));
        int last_extent = 0; //The extent which contains the last block
        int blocks = (int) (parent->size / A1FS_BLOCK_SIZE) + 1; //How many blocks this inode has allocated
        for (int i = 0; i < 512; i++) {
            blocks -= parent_extent_table->i_extents[i].count;
            if (blocks <= 0) {
                last_extent = i;
                break;
            }
        }
        int offset = (int) (((parent_extent_table->i_extents[last_extent].count - 1) * A1FS_BLOCK_SIZE) + space_used);
        int last_block = (int) (parent_extent_table->i_extents[last_extent].start +
                                parent_extent_table->i_extents[last_extent].count) - 1;
        printf("last extent start: %d \n", parent_extent_table->i_extents[last_extent].start);
        printf("offset: %d \n", offset);
        printf("last block: %d \n", last_block);

        //Now we can add the dir entry
        a1fs_dentry *dentry = img + (last_block * A1FS_BLOCK_SIZE) + offset;
        strncpy(dentry->name, file_name_actual, A1FS_NAME_MAX);
        dentry->ino = child_inode_number;

        //Now we update the parent's info
        parent->size = parent->size + sizeof(a1fs_dentry);
        clock_gettime(CLOCK_REALTIME, &parent->mtime);

        //Now we update the child's info
        child->size = 0;
        child->links = 2;
        clock_gettime(CLOCK_REALTIME, &child->mtime);
        child->mode = mode;

        flipBit(inodebitmap, child_inode_number);
        superblock->free_inodes_count -= 1;
    }

    (void)path;
    (void)mode;

    free(path_copy);
    return 0;
}


int rmdir2(const char* path, void * image)
{

	a1fs_superblock* superblock = (a1fs_superblock*)(image + A1FS_BLOCK_SIZE);
	a1fs_blk_t inode_table_block = superblock->inode_table;
	a1fs_inode* inodes = (a1fs_inode*)(image + (A1FS_BLOCK_SIZE * inode_table_block));
	a1fs_blk_t* blockbitmap = (a1fs_blk_t*)(image + (superblock->block_bitmap * A1FS_BLOCK_SIZE));
	a1fs_blk_t* inodebitmap = (a1fs_blk_t*)(image + (superblock->inode_bitmap * A1FS_BLOCK_SIZE));

	//Must first check if legal (namely the dir we want to remove is empty)
	int* path_info = path_lookup(path, image);
	int child_inode_number = path_info[1];
	a1fs_inode * child = &inodes[child_inode_number];
	if (child->size != 0) { return -1; }

	//Now we need to get the parent inode from which we want to remove the empty dir
	int parent_inode_number = path_info[3];
	a1fs_inode* parent = &inodes[parent_inode_number];
	a1fs_extent_table* parent_extent_table = (a1fs_extent_table*)(image + (parent->extent_table * A1FS_BLOCK_SIZE));
	printf("parent %d \n", parent_inode_number);

	// We want to know which index in parent to delete
	int child_index = path_info[4];
	printf("want to delete the %d dir from parent \n", child_index);

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

	printf("blocks passed %d \n", blocks_passed);
	int dest_offset = (sizeof(a1fs_dentry) * child_index)-(blocks_passed*A1FS_BLOCK_SIZE);
	printf("dest offset %d \n", dest_offset);

	//Scroll to very last extent/a1fs_dentry which we can use to "fill" the hole
	int src_offset = parent->size % A1FS_BLOCK_SIZE;
	src_offset += (parent_extent_table->i_extents[last_extent].count - 1) * A1FS_BLOCK_SIZE -sizeof(a1fs_dentry);
	printf("src offset %d \n", src_offset);
	printf("curr extent %d \n", curr_extent);
	printf("last extent %d \n", last_extent);
	if (last_extent == curr_extent && src_offset == dest_offset) { //We are removing the last dir
		parent->size -= sizeof(a1fs_dentry);
		printf("case 1\n");
		if (parent->size % A1FS_BLOCK_SIZE == 0) { //We shrunk the extent by one block
			printf("case 2 \n");
			parent_extent_table->i_extents[last_extent].count -= 1;
			flipBit(blockbitmap, parent_extent_table->i_extents[last_extent].start + parent_extent_table->i_extents[last_extent].count);
			superblock->free_blocks_count += 1;
		}
		if (parent->size == 0) {
			printf("case 4 \n");
			flipBit(blockbitmap, parent->extent_table);
			superblock->free_blocks_count += 1;
		}
		flipBit(inodebitmap, child_inode_number);
		superblock->free_inodes_count += 1;
		return 0;
	}
	printf("case 3 \n");
	memcpy(image + (parent_extent_table->i_extents[curr_extent].start * A1FS_BLOCK_SIZE) + dest_offset, image + (parent_extent_table->i_extents[last_extent].start * A1FS_BLOCK_SIZE) + src_offset, sizeof(a1fs_dentry));
	parent->size -= sizeof(a1fs_dentry);
	if (parent->size % A1FS_BLOCK_SIZE == 0) { //We shrunk the extent by one block
		printf("case 2 \n");
		parent_extent_table->i_extents[last_extent].count -= 1;
		flipBit(blockbitmap, parent_extent_table->i_extents[last_extent].start + parent_extent_table->i_extents[last_extent].count);
		superblock->free_blocks_count += 1;
	}
	flipBit(inodebitmap, child_inode_number);
	return 0;
}

int write2(const char *path, const char *buf, size_t size,
                      off_t offset, void * image)
{

    //TODO: write data from the buffer into the file at given offset, possibly
    // "zeroing out" the uninitialized range
    int* path_info = path_lookup(path, image);
    //a1fs_ino_t inode_num = (a1fs_ino_t) path_info[1];
    a1fs_superblock* superblock = (a1fs_superblock*)(image + A1FS_BLOCK_SIZE);
    a1fs_blk_t inode_table_block = superblock->inode_table;
    a1fs_inode* inodes = (a1fs_inode*)(image + (A1FS_BLOCK_SIZE * inode_table_block));
    a1fs_blk_t* blockbitmap = (a1fs_blk_t*)(image + (superblock->block_bitmap * A1FS_BLOCK_SIZE));
    //a1fs_blk_t* inodebitmap = (a1fs_blk_t*)(image + (superblock->inode_bitmap * A1FS_BLOCK_SIZE));
    a1fs_blk_t* data_blocks = (a1fs_blk_t*)(image + (superblock->data_blocks * A1FS_BLOCK_SIZE));
    a1fs_blk_t blocks_needed = 0;


    int inode_number = path_info[1];
    a1fs_inode* this_ino = &inodes[inode_number];
    a1fs_extent_table* extent_table = (a1fs_extent_table*)(image + (this_ino->extent_table * A1FS_BLOCK_SIZE));

    if ((inodes[path_info[1]].size == 0)) {
        printf("case 0 \n");
        superblock->free_blocks_count -= 1;
        if(superblock->free_blocks_count<0) {
            superblock->free_blocks_count += 1;
            return -1;
        }

        int i = findZeroBit(blockbitmap, superblock->blocks_count); //is non negative thx to previous check
        this_ino->extent_table = i;
        extent_table = (a1fs_extent_table*)(image + (this_ino->extent_table * A1FS_BLOCK_SIZE));
        flipBit(blockbitmap, i);

    }

    if ((inodes[path_info[1]].size == 0)) { // if the file size is 0
        char *words = create_string(buf, offset, 0);
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
            memcpy(image + (extent_info[0] * A1FS_BLOCK_SIZE), &words[copy_start], extent_count * A1FS_BLOCK_SIZE);
            printf("start: %d, size: %d \n", extent_info[0], extent_info[1]);
            printf("%d, %s, %d\n", extent_info[0], words, extent_count);
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
                memcpy(image + (extent_table->i_extents[ext_index].start * A1FS_BLOCK_SIZE) + (new_offset - curr_extent_size), &words[copy_index], extent_table->i_extents[ext_index].count * A1FS_BLOCK_SIZE);
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
                    memcpy(image + ((extent_table->i_extents[ext_index].start+extent_table->i_extents[last_ext].count-1) * A1FS_BLOCK_SIZE), &words[copy_index],  A1FS_BLOCK_SIZE);
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
                    memcpy(image + (extent_info[0] * A1FS_BLOCK_SIZE), &words[copy_index], extent_count * A1FS_BLOCK_SIZE);
                    copy_index += extent_count * A1FS_BLOCK_SIZE;
                    remaining_size -= extent_count*A1FS_BLOCK_SIZE;
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
    return 0;
}

int read2(const char *path, char *buf, size_t size, off_t offset, void * image)
{
    //TODO: read data from the file at given offset into the buffer
    a1fs_superblock* superblock = (a1fs_superblock*)(image + A1FS_BLOCK_SIZE);
    a1fs_blk_t inode_table_block = superblock->inode_table;
    a1fs_inode* inodes = (a1fs_inode*)(image + (A1FS_BLOCK_SIZE * inode_table_block));
    a1fs_blk_t* blockbitmap = (a1fs_blk_t*)(image + (superblock->block_bitmap * A1FS_BLOCK_SIZE));

    //Get file we want to truncate
    int* path_info = path_lookup(path, image);
    int child_inode_number = path_info[1];
    a1fs_inode* file = &inodes[child_inode_number];
    a1fs_extent_table* file_extent_table = (a1fs_extent_table*)(image + (file->extent_table * A1FS_BLOCK_SIZE));
    int last_ext = lastExtentIndex(file->size, file_extent_table);
    if(offset >= file->size){return 0;}
    int curr_size = 0;
    int start_found = 0;
    int amount_copied = 0;
    for(int i = 0; i <= last_ext; i++){
        file_extent_table->i_extents[i].count * A1FS_BLOCK_SIZE;
        printf("ext start: %d \n", file_extent_table->i_extents[i].start);
        if(start_found == 1){
            memcpy(&buf[amount_copied],file_extent_table->i_extents[i].start * A1FS_BLOCK_SIZE + image, file_extent_table->i_extents[i].count * A1FS_BLOCK_SIZE);
            amount_copied += (file_extent_table->i_extents[i].count * A1FS_BLOCK_SIZE);
        }
        if(offset >= curr_size && offset < curr_size+file_extent_table->i_extents[i].count * A1FS_BLOCK_SIZE){
            memcpy(buf,file_extent_table->i_extents[i].start * A1FS_BLOCK_SIZE + offset + image, file_extent_table->i_extents[i].count * A1FS_BLOCK_SIZE);
            printf("r: %d, %s, %d", file_extent_table->i_extents[i].start, buf, file_extent_table->i_extents[i].count * A1FS_BLOCK_SIZE);
            amount_copied += ((file_extent_table->i_extents[i].count * A1FS_BLOCK_SIZE)-(offset));
            start_found = 1;
        }
        curr_size += file_extent_table->i_extents[i].count * A1FS_BLOCK_SIZE;
    }
    return amount_copied;
}

//DEBUG END -----------------------------------------------------------------------------------------------------------

/**
 * Format the image into a1fs.
 *
 * NOTE: Must update mtime of the root directory.
 *
 * @param image  pointer to the start of the image.
 * @param size   image size in bytes.
 * @param opts   command line options.
 * @return       true on success;
 *               false on error, e.g. options are invalid for given image size.
 */
static bool mkfs(void* image, size_t size, mkfs_opts* opts)
{
	//TODO: initialize the superblock and create an empty root directory
	//NOTE: the mode of the root directory inode should be set to S_IFDIR | 0777

	a1fs_superblock* superblock = (a1fs_superblock*)(image + A1FS_BLOCK_SIZE);
	superblock->inodes_count = opts->n_inodes;
	superblock->free_inodes_count = opts->n_inodes;
	superblock->blocks_count = (size / A1FS_BLOCK_SIZE);
	int free_blocks = superblock->blocks_count - 2; //One for superblock and one for bios reserved block

	//Hypothetical partition of the disk
	int length_of_inode_table = (int)((opts->n_inodes * 64) / A1FS_BLOCK_SIZE) + 1;
	int length_of_inode_bitmap = (int)((opts->n_inodes) / A1FS_BLOCK_SIZE) + 1;

	free_blocks -= length_of_inode_table;
	free_blocks -= length_of_inode_bitmap;

	int length_of_block_bitmap = (int)(superblock->free_blocks_count / A1FS_BLOCK_SIZE) + 1;
	free_blocks -= length_of_block_bitmap;

	//If the partition is invalid we must return false
	if (free_blocks < 0) {
		return false;
	}
	superblock->free_blocks_count = free_blocks;

	//Real partition of the disk
	superblock->block_bitmap = 2;
	superblock->inode_bitmap = 2 + length_of_block_bitmap;
	superblock->inode_table = superblock->inode_bitmap + length_of_inode_bitmap;
	superblock->data_blocks = superblock->inode_table + length_of_inode_table;

	//Create root dir (Since its empty we wont allocate an extent table block)
	a1fs_inode* root_dir = (a1fs_inode*)(image + superblock->inode_table * A1FS_BLOCK_SIZE);
	root_dir->mode = (S_IFDIR | 0777);
	root_dir->links = 2;
	root_dir->size = 0;
	superblock->free_inodes_count -= 1;
	if (clock_gettime(CLOCK_REALTIME, &root_dir->mtime) != 0) {
		return false;
	}
	a1fs_blk_t * inodebitmap = (a1fs_blk_t *)(image + (superblock->inode_bitmap * A1FS_BLOCK_SIZE));
	flipBit(inodebitmap, 0);

	//We must also update the blocks bitmap
	a1fs_blk_t* blockbitmap = (a1fs_blk_t*)(image + (superblock->block_bitmap * A1FS_BLOCK_SIZE));
	flipBit(blockbitmap, 0); //For the bios reserved block
	flipBit(blockbitmap, 1); //For the superblock
	int i = 2;
	for (int k = 0; k < length_of_block_bitmap; k++) {
		flipBit(blockbitmap, i);
		i++;
	}
	for (int k = 0; k < length_of_inode_bitmap; k++) {
		flipBit(blockbitmap, i);
		i++;
	}
	for (int k = 0; k < length_of_inode_table; k++) {
		flipBit(blockbitmap, i);
		i++;
	}

	const char* path1 = "/dir1";
    const char* path2 = "/dir2";

    //printf("%d \n", blockbitmap[0]);
    //int* stuff = longest_extent(blockbitmap, superblock->blocks_count);
    //printf("%d, %d \n", stuff[0], stuff[1]);
    //char* str = create_string("kaan", 5, 0);
    //printf("%s \n", str);

	//const char* path3 = "/dir3";
	//const char* path4 = "/dir4";
	//const char* path5 = "/dir5";
	//const char* path6 = "/dir6";
	//const char* path7 = "/dir7";
	//const char* path8 = "/dir8";
	//const char* path9 = "/dir9";
	//const char* path10 = "/dir10";
	//const char* path11 = "/dir11";
	//const char* path12 = "/dir12";
	//const char* path13 = "/dir13";
	//const char* path14 = "/dir14";
	//const char* path15 = "/dir15";
	//const char* path16 = "/dir16";
	//const char* path17 = "/dir17";
	//const char* path18 = "/dir18";
	//mkdir2(path1, image);
	//mkdir2(path2, image);

	//rmdir2(path1, image);
	//mkdir2(path3, image);
	//mkdir2(path4, image);
	//mkdir2(path5, image);
	//mkdir2(path6, image);
	//mkdir2(path7, image);
	//mkdir2(path8, image);
	//mkdir2(path9, image);
	//mkdir2(path10, image);
	//mkdir2(path11, image);
	//mkdir2(path12, image);
	//mkdir2(path13, image);
	//mkdir2(path14, image);
	//mkdir2(path15, image);
	//mkdir2(path16, image);
	//mkdir2(path17, image);
	//mkdir2(path18, image);

    //printf("bitmap: %d \n", blockbitmap[0]);
	//create2(path1, S_IFREG | S_IRUSR, image);
    //printf("bitmap: %d \n", blockbitmap[0]);
    //write2(path1, "kaan", 4, 5, image);
    //write2(path1, "cinar", 5, 12, image);
    //printf("bitmap: %d \n", blockbitmap[0]);
    //char* str;
    //read2(path1, str, 4, 0, image);
    //printf("got: %s \n", str);


	superblock->magic = A1FS_MAGIC; //shows that the a1fs system was created on this disk.
	return true;

}


int main(int argc, char *argv[])
{
	mkfs_opts opts = {0};// defaults are all 0
	if (!parse_args(argc, argv, &opts)) {
		// Invalid arguments, print help to stderr
		print_help(stderr, argv[0]);
		return 1;
	}
	if (opts.help) {
		// Help requested, print it to stdout
		print_help(stdout, argv[0]);
		return 0;
	}

	// Map image file into memory
	size_t size;
	void *image = map_file(opts.img_path, A1FS_BLOCK_SIZE, &size);
	if (image == NULL) return 1;

	// Check if overwriting existing file system
	int ret = 1;
	if (!opts.force && a1fs_is_present(image)) {
		fprintf(stderr, "Image already contains a1fs; use -f to overwrite\n");
		goto end;
	}

	if (opts.zero) memset(image, 0, size);
	if (!mkfs(image, size, &opts)) {
		fprintf(stderr, "Failed to format the image\n");
		goto end;
	}

	ret = 0;
end:
	munmap(image, size);
	return ret;
}






