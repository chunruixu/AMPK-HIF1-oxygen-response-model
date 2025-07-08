getwd()    #查询当前目录
setwd("C:/Users/许春蕊/OneDrive - zzu.edu.cn/低氧投稿-补充在相应文件夹/GSE235562") #进入目标目录
getwd()    #查询当前目录

#section 1----################## calculate the average value of duplicate ##### (在此之前，删除无genename的行，在第一列加入id-1,2,3,...)

#load data
database_0 <- read.table(file = "0OriginalData_withdup.txt", 
                         sep = "\t", quote = "",
                         header = T, 
                         row.names = 1)

#calculate mean
database_mean <- aggregate(.~GENE_SYMBOL,mean,data = database_0)

#save file
#write.table(database_mean, file = "0database_mean.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


#section 2 ########################绘制差异显著箱线图###############
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

##总的数据
data <- read.table("1genelevel2.txt", header = T)#, row.names = 1, sep ="\t")
rownames(data) <- data[,1]
data <- data[,-1]


#单个基因的数据
data1 <- data[73:96,]  #SIRT -- 根据目标函数改变数字


#转factor
data1$sample <- factor(data1$sample, levels = c("Normoxia", "Hyperoxia"))

# 绘制箱线图
p <- ggplot(data1, aes(x = sample, y = value)) +
  geom_boxplot(aes(fill = sample), position = position_dodge(1),
               color = "black", size = 0.5) +
  # 加散点
  # geom_jitter(color = "#a6cee3", size = 2.5, alpha = 0.8) +
  stat_summary(
    fun = mean,
    geom = "point",
    aes(group = sample),
    position = position_dodge(1),
    color = "black",
    size = 2
  ) +
  stat_compare_means(          #标记会根据p值自动显示为ns（不显著）、*（p < 0.05）、**（p < 0.01）、***（p < 0.001）
    aes(x = sample),
    method = "t.test",
    label = "p.signif",
    label.y = max(data1$value) * 1.1,  # 调整标记的高度
    hjust = 1,  # 将标记向右对齐到 Hypoxia
    method.args = list(var.equal = TRUE)  # 指定等方差假设
  ) +
  scale_fill_manual(
    #values = c("Normoxia" = "lightblue", "Hypoxia" = "lightcoral"),
    values = c("Normoxia" = "white", "Hyperoxia" = "#FFA500")
  ) +
  #scale_x_discrete(limits = c("Normoxia", "Hypoxia")) +  # 确保Normoxia在左，Hypoxia在右
  labs(
    title = "GSE235562 Hypoxia data",
    x = "",
    y = "",
    fill = ""
    #x = "sample",
    #y = "expression level",
    #fill = "sample"
  ) +
  theme(
    panel.background = element_blank(),  # 去掉背景
    panel.grid = element_blank(),        # 去掉网格线
    plot.title = element_text(hjust = 0.5, size = 32),  # 增加标题字体大小
    axis.title = element_text(size = 30),  # 增加轴标题字体大小
    axis.text = element_text(size = 30),   # 增加轴刻度标签字体大小
    legend.text = element_text(size = 30), # 增加图例文本字体大小
    legend.position = "top",
    axis.line = element_line(),          # 添加坐标轴线
    axis.ticks = element_line()          # 添加坐标轴刻度线
  ) +
  guides(fill = guide_legend(keywidth = 3, keyheight = 3)) # 增加图例符号的大小


# 显示图形
print(p)

# 保存为高分辨率PNG文件
ggsave(
  filename = "D:/我的坚果云/俄罗斯高低氧/2024ACM - 论文1/投稿/SIRT1 level hyperoxia.png",  # 指定保存路径和文件名
  plot = p,
  width = 8,  # 图片宽度（英寸）
  height = 8,  # 图片高度（英寸）
  dpi = 300,   # 分辨率（每英寸点数）
  units = "in" # 单位（英寸）
)


#section 3-----#################火山图#########s##########
#load data - EXCEL OR resdata from DEGseq2 --需加入COL or Change 列
resdata <- read.table(file = "1UPandDown_excel_controlVSHOX80_p0.05_FC1.1and0.9.txt", 
                      sep = "\t", 
                      header = T, 
                      row.names = 1)

#删掉#DIV/0!



library(ggplot2)
library(dplyr)
library(ggrepel)
# 使用geom_label_repel函数为火山图添加目标基因的标签
# 读入目标基因列表
gene_label <- read.table("1genelist1.txt",header = T)
# 转换为数据框格式
gene_label <- data.frame(gene_label)
# 为数据框添加新的列，等于目标基因列表，方便下一步left_join合并数据
gene_label$list <- gene_label$label
# dplyr包中的left_join函数合并resdata数据表格和gene_label表格
#resdata2 <- resdata %>% left_join(gene_label, by=c("Row.names"="label"))
resdata2 <- resdata %>% left_join(gene_label, by=c("GENE_SYMBOL"="label"))#在最后加一列label和GENE_SYMBOL对应

#打开已存文档
#dat = read.delim(file.choose(), head = T, row.names = 1)

pdf("VolcanoPlot1.pdf",width = 7,height = 5)
p <- ggplot(data=resdata2,                 #绘图所使用的数据
            mapping=aes(x=logFC, #X轴使用的数据
                        y=P,   #Y轴使用的数据
                        color=COL))+   #数据呈现颜色按change列区分（共3种）
  geom_point(size=2)+                      #绘制散点图 
  scale_color_manual(values =c("#546de5", "#d2dae2","#ff4757"))+    #设定3种散点的颜色
  #xlim(-5,5)+ylim(0,80)+        #X轴和Y轴上下限
  #xlim(min(resdata2$logFC), max(resdata2$logFC)) +
  xlim(-max(resdata2$logFC), max(resdata2$logFC)) +
  #xlim(min(resdata2$logFC), -min(resdata2$logFC)) +
  ylim(min(resdata2$P), max(resdata2$P))+
  labs(x="log2(fold change)",       #命名X轴标题
       y="-log10(p-value)")+          #命名Y轴标题
  geom_hline(yintercept=-log10(0.05),linetype=4)+   #按照筛选差异表达基因的标准padj<0.05添加虚线（虚线线型：4）
  geom_vline(xintercept=c(-1,1),linetype=4)+        #按照筛选差异表达基因的标准|log2FC|>1添加虚线
  theme_bw()+  #ggplot2默认theme_grey主题（灰白主题），使用theme_bw()变为黑白主题
  theme(axis.text=element_text(size=20),     #设置坐标轴字体大小
        axis.title=element_text(size=20),    #设置坐标轴标题字体大小
        legend.text=element_text(size=15),   #设置图例字体大小
        legend.title=element_text(size=15))+ #设置图例标题字体大小
  geom_label_repel(aes(label = list),
                   size = 6,#3,
                   #box.padding = unit(0.5, "lines"),point.padding = unit(0.5, "lines"), #设置散点与标签文本框间距
                   segment.color = "black",  #连接点与标签的线段的颜色
                   show.legend =F)  #避免图例显示问题

print(p)


# 保存为高分辨率PNG文件
ggsave(
  filename = "D:/我的坚果云/俄罗斯高低氧/2024ACM - 论文1/投稿/volcano_hypoxia2.png",  # 指定保存路径和文件名
  plot = p,
  width = 8,  # 图片宽度（英寸）
  height = 8,  # 图片高度（英寸）
  dpi = 300,   # 分辨率（每英寸点数）
  units = "in" # 单位（英寸）
)

dev.off()



#section 4######################GO富集分析的可视化###################
##富集分析表格用DAVID和g profiler做完
library (dplyr)
library (ggplot2)  
#install.packages("tidyverse")
library(tidyverse)
#install.packages("openxlsx")
library(openxlsx)


#数值导入#
BP_UP = read.xlsx('2富集分析.xlsx',sheet= "BP_UP_selected",sep=',')
CC_UP = read.xlsx('2富集分析.xlsx',sheet= "CC_UP_selected",sep=',')
MF_UP = read.xlsx('2富集分析.xlsx',sheet= "MF_UP_selected",sep=',')
BP_ALL = read.xlsx('2富集分析.xlsx',sheet= "GOBP_TOP100",sep=',')
CC_ALL = read.xlsx('2富集分析.xlsx',sheet= "GOCC_TOP100",sep=',')
MF_ALL = read.xlsx('2富集分析.xlsx',sheet= "GOMF_TOP100",sep=',')

#head(BP_UP)
BP = BP_UP
CC = CC_UP
MF = MF_UP
BP = BP_ALL
CC = CC_ALL
MF = MF_ALL
BP = separate(BP,Term, sep="~",into=c("ID","Description"))
CC = separate(CC,Term, sep="~",into=c("ID","Description"))
MF = separate(MF,Term, sep="~",into=c("ID","Description"))
colnames(BP)[5] <- "Gene_Ratio"
colnames(CC)[5] <- "Gene_Ratio"
colnames(MF)[5] <- "Gene_Ratio"

#提取各组数据需要展示的数量#
display_number = c(4, 3, 4)  ##这三个数字分别代表选取的BP、CC、MF的数量###需要改!!!!!!!!!!!
go_result_BP = as.data.frame(BP)[1:display_number[1], ]
go_result_CC = as.data.frame(CC)[1:display_number[2], ]
go_result_MF = as.data.frame(MF)[1:display_number[3], ]


#将提取的各组数据进行整合
go_enrich = data.frame(
  ID=c(go_result_BP$ID, go_result_CC$ID, go_result_MF$ID),  #指定ego_result_BP、ego_result_CC、ego_result_MFID为ID                        
  Description=c(go_result_BP$Description,go_result_CC$Description,go_result_MF$Description),
  GeneNumber=c(go_result_BP$Count, go_result_CC$Count, go_result_MF$Count), #指定ego_result_BP、ego_result_CC、ego_result_MF的Count为GeneNumber
  P=c(go_result_BP$P,go_result_CC$P,go_result_MF$P),
  type=factor(c(rep("Biological Process", display_number[1]), #设置biological process、cellular component、molecular function 的展示顺序
                rep("Cellular Component", display_number[2]),
                rep("Molecular Function", display_number[3])),
              levels=c("Biological Process", "Cellular Component","Molecular Function" )))



##设置GO term名字长短，过长则设置相应长度
for(i in 1:nrow(go_enrich)){
  description_splite=strsplit(go_enrich$Description[i],split = " ")
  description_collapse=paste(description_splite[[1]][1:5],collapse = " ")  #选择前五个单词GO term名字
  go_enrich$Description[i]=description_collapse
  go_enrich$Description=gsub(pattern = "NA","",go_enrich$Description)  #gsub()可以用于字段的删减、增补、替换和切割，可以处理一个字段也可以处理由字段组成的向量。gsub(“目标字符”, “替换字符”, 对象)
}



#转成因子，防止重新排列
#go_enrich$type_order=factor(rev(as.integer(rownames(go_enrich))),labels=rev(go_enrich$Description))
go_enrich$type_order = factor(go_enrich$Description,levels=go_enrich$Description,ordered = T)

#head(go_enrich)

#纵向柱状图#
p <- ggplot(go_enrich,
            aes(x=type_order,y=P, fill=type)) + #x、y轴定义；根据type填充颜色
  geom_bar(stat="identity", width=0.8) +  #柱状图宽度
  scale_fill_manual(values = c("#6666FF", "#33CC33", "#FF6666") ) +  #柱状图填充颜色
  coord_flip() +  #让柱状图变为纵向
  xlab("GO term") +  #x轴标签
  ylab("-log10(p-value)") +  #y轴标签
  labs(title = "GO Terms Enrich")+  #设置标题
  theme_bw()

print(p)


#纵向柱状图#
go_enrichBPandMF <- go_enrich[c(1:5, 11:15), ]#HYPOXIA
go_enrichBPandMF <- go_enrich[c(1:4, 8:11), ]#HYPEROXIA

# 按 type 分组，然后在每个组内按 P 值从大到小排序
#go_enrichBPandMF <- go_enrichBPandMF %>%
#  arrange(type, desc(P)) %>%
#  mutate(type_order = factor(type_order, levels = type_order))


p <- ggplot(go_enrichBPandMF,
            aes(x=type_order,y=P, fill=type)) + #x、y轴定义；根据type填充颜色
  geom_bar(stat="identity", width=0.8) +  #柱状图宽度
  scale_fill_manual(values = c("#1E73FF", "#33CC33") ) +  #柱状图填充颜色
  coord_flip() +  #让柱状图变为纵向
  xlab("GO term") +  #x轴标签
  ylab("-log10(p-value)") +  #y轴标签
  labs(title = "GO Terms Enrich")+  #设置标题
  theme_bw()



# 按 type 分组，然后在每个组内按 P 值从小到大排序
go_enrichBPandMF <- go_enrichBPandMF %>%
  arrange(type, P) %>%
  mutate(type_order = factor(type_order, levels = type_order))

# 绘制横向柱状图
p <- ggplot(go_enrichBPandMF, aes(x = type_order, y = P, fill = type)) +
  geom_bar(stat = "identity", width = 0.8) +  # 柱状图宽度
  scale_fill_manual(values = c("#1E73FF", "#33CC33")) +  # 柱状图填充颜色
  coord_flip() +  # 翻转坐标轴，使柱状图变为横向
  xlab("GO term") +  # x轴标签
  ylab("-log10(p-value)") +  # y轴标签
  labs(title = "GO Terms Enrich") +  # 设置标题
  theme(axis.title.x = element_text(face = "bold",size = 16),
        axis.title.y = element_text(face = "bold",size = 16),
        legend.title = element_text(face = "bold",size = 16))+
  theme_bw()

print(p)



# 保存为高分辨率PNG文件
ggsave(
  filename = "D:/我的坚果云/俄罗斯高低氧/2024ACM - 论文1/投稿/GO_hyperoxia.png",  # 指定保存路径和文件名
  plot = p,
  width = 8,  # 图片宽度（英寸）
  height = 8,  # 图片高度（英寸）
  dpi = 300,   # 分辨率（每英寸点数）
  units = "in" # 单位（英寸）
)



#横向柱状图#
pdf("GO_bar.pdf",width = 10,height = 5)
ggplot(go_enrich,
       aes(x=type_order,y=GeneNumber, fill=type)) +  #x、y轴定义；根据type填充颜色
  geom_bar(stat="identity", width=0.8) +  #柱状图宽度
  scale_fill_manual(values = c("#6666FF", "#33CC33", "#FF6666") ) + #柱状图填充颜色
  xlab("GO term") + #x轴标签
  ylab("Gene_Number") +  #y轴标签
  labs(title = "GO Terms Enrich")+ #设置标题
  theme_bw() +
  theme(axis.text.x=element_text(family="sans",face = "bold", color="gray50",angle = 70,vjust = 1, hjust = 1 )) #对字体样式、颜色、还有横坐标角度（）
dev.off()


#分组柱状图绘制# ###分组：选BP or CC or MF
#pdf("BPALL_bar.pdf",width = 10,height = 5)
target <- go_result_BP
#PValue转成数值类型
target$PValue <- as.numeric(target$PValue)
ggplot(target,aes(x = Gene_Ratio,y = Description,fill=PValue))  +  #fill=PValue，根据PValue填充颜色；
  geom_bar(stat="identity",width=0.8 ) + ylim(rev(go_result_BP$Description)) + #柱状图宽度与y轴显示顺序
  theme_test() + #调整背景色
  scale_fill_gradient(low="red",high ="blue") +  #调节填充的颜色
  labs(size="Genecounts",x="GeneRatio",y="GO term",title="GO_BP") + #设置图内标签（x、y轴、标题）
  theme(axis.text=element_text(size=10,color="black"),
        axis.title = element_text(size=16),title = element_text(size=13))
dev.off()

#气泡图#
ego = rbind(go_result_BP,go_result_CC,go_result_MF)  
ego = as.data.frame(ego)
rownames(ego) = 1:nrow(ego)
ego$order=factor(rev(as.integer(rownames(ego))),labels = rev(ego$Description))

head(ego)

#转数据类型
ego$PValue <- as.numeric(ego$PValue)

ggplot(ego,aes(y=order,x=Gene_Ratio))+
  geom_point(aes(size=Count,color=PValue))+
  scale_color_gradient(low = "red",high ="blue")+
  labs(color=expression(PValue,size="Count"),
       x="Gene Ratio",y="GO term",title="GO Enrichment")+
  theme_bw()



#分组气泡图绘制#
pdf("BPALL_bubble.pdf",width = 10,height = 5)
ggplot(go_result_BP,  #根据所绘制的分组进行绘制
       aes(x=Gene_Ratio,y=Description,color= PValue)) +  #color= -1*PValue可定义成倒数
  geom_point(aes(size=Count))  +
  ylim(rev(go_result_BP$Description)) +  #rev相反的意思
  labs(size="Genecounts",x="GeneRatio",y="GO_term",title="GO_BP") +  #调标签名称
  scale_color_gradient(low="red",high ="blue") +   #改颜色
  theme(axis.text=element_text(size=10,color="black"),  
        axis.title = element_text(size=16),title = element_text(size=13))+
  theme_bw()
dev.off()

pdf("MFUP_bubble.pdf",width = 7,height = 5)
ggplot(go_result_MF,  #根据所绘制的分组进行绘制
       aes(x=Gene_Ratio,y=Description,color= PValue)) +  #color= -1*PValue可定义成倒数
  geom_point(aes(size=Count))  +
  ylim(rev(go_result_MF$Description)) +  #rev相反的意思
  labs(size="Genecounts",x="GeneRatio",y="GO_term",title="GO_MF") +  #调标签名称
  scale_color_gradient(low="red",high ="blue") +   #改颜色
  theme(axis.text=element_text(size=10,color="black"),  
        axis.title = element_text(size=16),title = element_text(size=13))+
  theme_bw()
dev.off()


######################KEGG富集分析的可视化#######################
#数据导入#
kk_result= read.xlsx('2富集分析.xlsx', sheet= "KEGG_TOP100", sep = ',')
kk_result = separate(kk_result,Term, sep=":",into=c("ID","Description"))
colnames(kk_result)[5] <- "Gene_Ratio"
#数据处理#
display_number = 6#显示数量设置
kk_result = as.data.frame(kk_result)[1:display_number[1], ]
kk = as.data.frame(kk_result)
rownames(kk) = 1:nrow(kk)
kk$order=factor(rev(as.integer(rownames(kk))),labels = rev(kk$Description))

#柱状图#
pdf("KEGGALL_bar.pdf",width = 7,height = 5)
p <- ggplot(kk,aes(y=order,x=P,fill=as.numeric(Count)))+
  geom_bar(stat = "identity",width=0.8)+ #柱状图宽度设置
  scale_fill_gradient(low = "#FFB6C1",high ="#FF0000" )+
  labs(title = "KEGG Pathways Enrichment",  #设置标题、x轴和Y轴名称
       x = "-log10(p-value)",
       y = "Pathway")+
  theme(axis.title.x = element_text(face = "bold",size = 16),
        axis.title.y = element_text(face = "bold",size = 16),
        legend.title = element_text(face = "bold",size = 16))+
  theme_bw()

print(p)
# 保存为高分辨率PNG文件
ggsave(
  filename = "D:/我的坚果云/俄罗斯高低氧/2024ACM - 论文1/投稿/GO_hypoxia2.png",  # 指定保存路径和文件名
  plot = p,
  width = 8,  # 图片宽度（英寸）
  height = 8,  # 图片高度（英寸）
  dpi = 300,   # 分辨率（每英寸点数）
  units = "in" # 单位（英寸）
)

dev.off()


#气泡图#
pdf("KEGGALL_bubble.pdf",width = 7,height = 5)
ggplot(kk,aes(y=order,x=Gene_Ratio))+
  geom_point(aes(size=Count,color=PValue))+
  scale_color_gradient(low = "red",high ="blue")+
  labs(color=expression(PValue,size="Count"),
       x="Gene Ratio",y="Pathways",title="KEGG Pathway Enrichment")+
  theme_bw()

p <- ggplot(kk, aes(y = reorder(Description, P), x = P, size = Count, color = P)) +
  geom_point(alpha = 0.7) +
  scale_size_continuous(range = c(5, 15)) +  # 设置点的大小范围
  scale_color_gradient(low = "#FFB6C1", high = "#FF0000") +  # 设置颜色渐变
  scale_x_continuous(limits = c(min(kk$P)*1.2, max(kk$P)*1.2)) +  # 设置x轴的范围
  labs(title = "KEGG Pathway Dot Plot", x = "-log10(p-value)", y = "Pathway",
       size = "Gene Count", color = "-log10(p-value)") +
  scale_x_continuous(limits = c(1, max(kk$P)*1.2)) +  # 设置x轴的范围从0开始
  theme(axis.title.x = element_text(face = "bold",size = 16),
        axis.title.y = element_text(face = "bold",size = 16),
        legend.title = element_text(face = "bold",size = 16))+
  theme_bw()

# p <- ggplot(kk, aes(y = Description, x = P, size = Count, color = P)) +
#   geom_point(alpha = 0.7) +
#   scale_size_continuous(range = c(5, 15)) +  # 设置点的大小范围
#   scale_color_gradient(low = "#FFB6C1", high = "#FF0000") +  # 设置颜色渐变
#   scale_x_continuous(limits = c(min(kk$P)*1.2, max(kk$P)*1.2)) +  # 设置x轴的范围
#   labs(title = "KEGG Pathway Dot Plot", x = "-log10(p-value)", y = "Pathway",
#        size = "Gene Count", color = "-log10(p-value)") +
#   theme(axis.title.x = element_text(face = "bold",size = 16),
#         axis.title.y = element_text(face = "bold",size = 16),
#         legend.title = element_text(face = "bold",size = 16))+
#   theme_bw()

print(p)

# 保存为高分辨率PNG文件
ggsave(
  filename = "D:/我的坚果云/俄罗斯高低氧/2024ACM - 论文1/投稿/KEGG_hyperoxia.png",  # 指定保存路径和文件名
  plot = p,
  width = 8,  # 图片宽度（英寸）
  height = 8,  # 图片高度（英寸）
  dpi = 300,   # 分辨率（每英寸点数）
  units = "in" # 单位（英寸）
)


dev.off()


#section 5#################3相关性分析################
# 安装并加载必要的包
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!requireNamespace("ggpubr", quietly = TRUE)) {
  install.packages("ggpubr")
}

library(ggplot2)
library(ggpubr)
library(openxlsx)

######皮尔森相关系数#####
#"data<-read.delim(file.choose(),row.names=1)"

data = read.xlsx('3关联分析.xlsx',sheet= "HIF1目的基因整理",sep=',')
# data <- na.omit(data)
# rownames(data) <- data$GENE_SYMBOL
# data <- data[, -1]
# data<-as.matrix(data)
# 
# length1<-length(data[1,])
# length2<-length(data[,1])
# x <- data[1,1:length1]
# x <- as.numeric(x)
# 
# 
# pvalue<-c()
# cor<-c()
# for (i in 2:length2)
# {
#   y<-data[i,1:length1]
#   y<-as.numeric(y)
#   p<-cor.test(x,y)$p.value
#   pvalue<-c(pvalue,p)
#   co<-cor(x,y)
#   cor<-c(cor,co)
# }
# pvalue<-c("NA",pvalue)
# cor<-c("NA",cor)
# correlation<-cbind(cor,pvalue)
# row.names(correlation)<-row.names(data)
# write.csv(correlation,"Pearson-cor-hif1a.csv")


#######斯皮尔曼相关系数########
data = read.xlsx('3关联分析.xlsx',sheet= "目的基因整理",sep=',')
data <- na.omit(data)
rownames(data) <- data$GENE_SYMBOL
data <- data[, -1]
data<-as.matrix(data)

length1<-length(data[1,])
length2<-length(data[,1])
x<-data[1,1:length1]
x<-as.numeric(x)
pvalue<-c()
cor<-c()
for (i in 2:length2)
{
  y<-data[i,1:length1]
  y<-as.numeric(y)
  p<-cor.test(x,y, method = 'spearman')$p.value
  pvalue<-c(pvalue,p)
  co<-cor(x,y, method = 'spearman')
  cor<-c(cor,co)
}
pvalue<-c("NA",pvalue)
cor<-c("NA",cor)
correlation<-cbind(cor,pvalue)
row.names(correlation)<-row.names(data)
write.csv(correlation,"spearman-cor-HIF1A.csv")