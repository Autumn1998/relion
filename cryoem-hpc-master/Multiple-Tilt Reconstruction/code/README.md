代码根目录：

41服务器 /home/wangzihao/CudaFBP/

建议首先看明白 multirecon   
main函数文件：main-fbp.cpp  
重构函数文件：atom-fbp.cpp  
其他文件是处理参数或者处理mrc读入输出的文件
 
以multirecon为例，代码有两种运行方式

1.使用nsight打开multirecon这个项目，然后在debug下面编译出可执行程序  
2.使用makefile文件，将代码编译出可执行程序

multirecon 使用最简单的BPT重构，即back-projection
算法参考链接:  
https://www.sciencedirect.com/topics/computer-science/backprojection

具体执行方式在multirecon目录下有sh脚本  
/home/wangzihao/CudaFBP/multirecon/singleexp/runmpi.sh

全部代码的运行方式在以下路径：   
/home/wangzihao/ysingleexp
/home/wangzihao/zdoubleexp