import os

def get_parser_file_type(parser, must_exist = False):

    def _file_type(path):
    
        path = os.path.expanduser(path)
    
        if must_exist:
            if not os.path.exists(path):
                parser.error('File doesn\'t exist: %s' % path)
            elif not os.path.isfile(path):
                parser.error('Not a file: %s' % path)
            else:
                return path
        else:
        
            dir_path = os.path.dirname(path)
        
            if not os.path.exists(dir_path):
                parser.error('Parent directory doesn\'t exist: %s' % dir_path)
            else:
                return path
    
    return _file_type

def get_parser_directory_type(parser, create_if_not_exists = False):
    
    def _directory_type(path):
    
        path = os.path.expanduser(path)
    
        if not os.path.exists(path):
            if create_if_not_exists:
                if os.path.exists(os.path.dirname(path)):
                    os.mkdir(path)
                    return path
                else:
                    parser.error('Cannot create empty directory (parent directory doesn\'t exist): %s' % path)
            else:
                parser.error('Path doesn\'t exist: %s' % path)
        elif not os.path.isdir(path):
            parser.error('Not a directory: %s' % path)
        else:
            return path
        
    return _directory_type