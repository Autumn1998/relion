+ 全部数据路径

41服务器 /home/wangzihao/UCSD_data

+ 双轴数据路径

/home/wangzihao/eel-1024

+ 以双轴为例，各个文件的简单介绍
1. eel-tomo4a.st  
mrc格式文件，是121张1024*1024的投影图片栈(所以是st结尾(stack))，a是指绕着一个轴进行旋转
2. eel-tomo4a.txbr  
文本文件，363行，每三行对应eel-tomo4a.st中的一个投影角度
3. eel-tomo4b.st  
b是指绕着与a垂直90度的方向旋转
4. eel-tomo4b.txbr   
与a同理 
5. eel-tomo4.st  
没有后缀的这个文件是将a与b的图片合并到一起，具有242张图片
6. eel-tomo4.txbr
eel-tomo4a.txbr 与eel-tomo4b.txbr 拼接而成的文本文件

st文件可以使用imod的两个指令查看  
指令一：
imod eel-tomo4a.st  
指令二：
header eel-tomo4a.st

补充内容：
imod的下载路径  
http://bio3d.colorado.edu/imod/

MRC图像格式介绍  
https://www.ccpem.ac.uk/mrc_format/mrc2014.php

