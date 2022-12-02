# Method 1, Google Driver

Sharing link: https://drive.google.com/drive/folders/1mgqjO5ABV5mQxYRUbZEp4XVwCJ6BKW61?usp=sharing

seu_A124.counts_anno.rds is the single cell RNA-seq data from 124 patients with liver cancer. 
seu_mouse_liverCancer.anno.rds is mouse liver cancer data. 

To avoid network interruption, we have splited the seurat objects into several 1 Gb large files. You can download them and combine them to generate the complete Rdata object. 
e.g.
```
cat seu_A124.counts_anno.rds.a* > seu_A124.counts_anno.rds
```

# Method 2, ftp download (recommendation for users in China)

IP: 123.124.148.211
Name: ftpguest1
Password: 123456
Recommended software: Filezilla
