# Navigating your directories on command line
A basic requirement to work on the command line is to be able to effectively move around your directory, create/delete/move directories and files, set permissions for directories and files.

## Where are you?
First off where are you? Are you in your home directory or another directory? To find your home directory path or to which directory you are currently in, there are some command line commands you should remember
```
# Your home directory path
echo $HOME

# The current directory you are in
pwd

# To move up one directory above
$ cd ../

# To move to home directory
$ cd ~/
```

## Creating, deleting, or moving directories and files
```
# To create a directory
$ mkdir #NAME_DIRECTORY

# Removing directory or files
$ rm -r #DIRECTORYNAME
$ rm #FILENAME
```

## Permissions
Permissions allows you or others to define who gets access to read, edit, or execute certain files.
```
# To see permission details, use ls with parameter -l (long list)
$ ls -l

-rw-r--r-- 4 root root 183 Month Day HH:MM filename.txt
```
