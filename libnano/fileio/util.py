import os
import mmap
 
 
def getFileContents(fp, max_in_mem_size=200):
    ''' Return the contents of file at `fp` either as a string (in memory) or
    as a memory-mapped file (if the file size is > `max_in_mem_size` in MB)
    '''
    contents_object = None
    is_mmap = False
    with open(fp, 'r+b') as fd:  # r+ read mode necessary for mmap
        if os.path.getsize(fp)/1024 > max_in_mem_size:
            contents_object = mmap.mmap(fd.fileno(), 0)  # 0->mmap entire file
            is_mmap = True
        else:
            contents_object = fd.read()
    return contents_object, is_mmap