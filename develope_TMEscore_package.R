


library(devtools)
library(roxygen2)

################################
#' @添加依赖包
use_package("dplyr")
use_package("ggplot2")
use_package("survival")
use_package("survminer")
use_package("tidyverse")
###############################



#' @把生成好的数据加载到最外部的工作环境下
##############################
(load("eset_acrg.RData"))
devtools::use_data(eset)
#############################



#' @常用快捷方式
#' @RStudio的快捷键来实现文档的编辑：Ctrl+Shift+Alt+R（光标放在函数名上--先写好function）


#' @把写好的function加载到程序中:Ctrl-Shift-L
devtools::load_all() # 把包骨架文件夹中的 R 文件夹中的所有 .R 文件读进来


#' @测试function
fbm() # 测试自己写的程序
fbm(hurst=0.2, n=1000) # 再测试自己写的程序
#########################


#' @生成对应的*.Rd文件在man文件夹中:Ctrl+Shift+D
devtools::document()


#' @打包
#' @就会在与包文件夹平行的文件夹中生成
#' somebm_0.1.tar.gz 类似的打包文件。可以在 R 环境中使用
#' install.packages('~/somebm_0.1.tar.gz', type='source') 来安装！
devtools:: build()


#' @DOCUMENTING YOUR PACKAGE
devtools::use_vignette('introduction')
#This will create a vignette/ folder and an R Markdown file, where you can mix text, code, and code output for a polished tutorial on using your package. You can access this vignette when the package is installed by using:
vignette('introduction', package = 'TMEscore')


# Unit testing (via the testthat package and
devtools::use_testthat() #allows you to include tests to confirm that your functions work as intended
#run Ctrl-Shift-T to run
devtools::test() #and run the unit tests).


#' @设置账号和账户
git config --global user.name "DongqiangZeng0808"
git config --global user.email "708902023@qq.com"


#' @参考链接
# https://pjnewcombe.wordpress.com/2014/05/31/using-rstudio-to-create-and-maintain-an-r-package-on-github/
#' @or create a new repository on the command line
echo "# TMEscore" >> README.md
git init
git add README.md
git commit -m "first commit"
git remote add origin git@github.com:DongqiangZeng0808/TMEscore.git
git push -u origin master

#' @or push an existing repository from the command line
git remote add origin git@github.com:DongqiangZeng0808/TMEscore.git
git push -u origin master
############################################



