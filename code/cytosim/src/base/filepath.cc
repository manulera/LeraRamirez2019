// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "filepath.h"
#include <sys/param.h>
#include <sys/stat.h>
#include <cstdlib>
#include <libgen.h>
#include <unistd.h>
#include <dirent.h>
#include "exceptions.h"


std::string FilePath::get_dir()
{
    char cwd[1024];
    getcwd(cwd, sizeof(cwd));
    return cwd;
}


bool FilePath::is_dir(std::string const& path)
{
    struct stat s;
    if ( 0 == stat(path.c_str(), &s) )
        return S_ISDIR(s.st_mode);
    return false;
}


int FilePath::change_dir(std::string const& path)
{
    if ( path == "." )
        return 0;
    return chdir(path.c_str());
}


std::vector<std::string> FilePath::list_dir(std::string const& path)
{
    std::vector<std::string> res;
    DIR * dp = opendir(path.c_str());
    if ( dp )
    {
        struct dirent * ep = readdir(dp);
        while ( ep )
        {
            res.push_back(ep->d_name);
            ep = readdir(dp);
        }
        
        closedir(dp);
    }
    return res;
}


std::vector<std::string>  FilePath::list_dir(std::string const& path, std::string const& ext)
{
    std::vector<std::string> res;
    DIR * dp = opendir(path.c_str());
    if ( dp )
    {
        struct dirent * ep = readdir(dp);
        while ( ep )
        {
            std::string name(ep->d_name);
            if ( name.size() > ext.size()
                &&  0==name.compare(name.size()-ext.size(), ext.size(), ext) )
                res.push_back(ep->d_name);
            ep = readdir(dp);
        }
        
        closedir(dp);
    }
    
#if ( 0 )
    for ( unsigned i = 0; i < res.size(); ++i )
        std::clog << "   " << res[i] << std::endl;
#endif
    return res;
}


std::string FilePath::dir_part(std::string const& path)
{
    char* res, tmp[PATH_MAX];
    
    if ( 0 == realpath(path.c_str(), tmp) )
        return ".";
    
    res = dirname(tmp);
    
    if ( 0 == res )
        throw InvalidIO("FilePath: stdlib::dirname() failed");
    
    return res;
}


bool FilePath::is_file(std::string const& path)
{
    struct stat s;
    if ( 0 == stat(path.c_str(), &s) )
        return S_ISREG(s.st_mode);
    return false;
}


std::string FilePath::file_part(std::string const& path)
{
    char * res = basename(const_cast<char*>(path.c_str()));
    
    if ( 0 == res )
        throw InvalidIO("FilePath: stdlib::basename() failed");
    
    return res;
}


std::string FilePath::full_name(std::string const& dir, std::string const& file)
{
    //if a full path is already specified, we do nothing
    if ( dir.size() > 0  &&  file.size() > 0  &&  file[0] != '/' )
    {
        std::string res = dir + file;
        
        //remove trailling '/' if present
        const int x = res.size()-1;
        if ( 0 <= x  &&  res[x] == '/' )
            res[x] = '\0';
        
        return res;
    }
    return file;
}

