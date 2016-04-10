# Usage

## 1. Downloading Massbank
```     
svn checkout http://www.massbank.jp/SVN/OpenData/record/
```
## 2. Copying all the individual Massbank records into allMBFiles folder
For windows :     
http://stackoverflow.com/questions/585091/using-xcopy-to-copy-files-from-several-directories-to-one-directory

For Mac/Unix :     
http://unix.stackexchange.com/questions/67503/move-all-files-with-a-certain-extension-from-multiple-subdirectories-into-one-di
## 3. use python script to generate mzML
```     
python parseMBFiles.py allMBFiles/
```

Note : For updating MB2HMDBmapping.csv you need to update HMDB first : 
https://github.com/epoyraz/updateHMDB
