# 基于相似度的*Needleman*-*Wunsch*序列对齐算法
### 0 背景引入

*Needleman-Wunsch*序列对齐算法最初源自于计算生物学，通过最大化局部相似性来对齐DNA序列、蛋白质序列等生物序列。比如输入两序列`ATCGCGGA`和`AGCTCAAT`:

<img src="sequence-alignment-algorithm-based-on-similarity-of-characters/DNA序列配对.png" alt="DNA序列配对"  />

同样在中文纠错、变体字识别等任务场景中也会出现需要使用序列对齐算法来快速完成序列标注，为接下来的模型训练做好准备工作。

![变体字序列对齐示例](sequence-alignment-algorithm-based-on-similarity-of-characters/变体字序列对齐示例.png)

传统的序列对齐算法的评价指标为对齐后两序列的最大公共字串长度最长(可不连续)或者两序列的编辑距离最短。但是由人为精心构造的变体序列是专门为了绕过审查监管的，在极端情况下，甚至该序列中每个字都与标注序列不同。在这种情况下，我们无法直接得到最长公共字符串，从而传统的序列对齐算法无能为力。

**因此，我在原有算法的基础上，提出了基于字符相似度的序列对齐算法，使得对齐后两序列各位置字符相似度之和最大。**



### 1 字符相似度计算

在序列中，原文本中字符应包括汉字、英文字母、数字、常见中英文标点符号以及unicode中其他非常字符。字符相似度计算函数$simi\_score(i, j)$用来给出两字符的相似度得分, 对于不同类型的字符间相似度计算规则如下：

1. 若两字符均为汉字，则采用基于汉字音形码的相似度计算函数，详见[这个博客](https://blog.csdn.net/chndata/article/details/41114771)。音形码可以表示汉字在读音上和字形上的特征，包括了声母、韵母、音调、汉字结构、四角编码和笔画数。将汉字使用音形码嵌入后通过特定的计算公式计算相似度。
2. 若两字符均为字母或者均为数字，则直接判断两字符是否相等。
3. 若两字符为中英文常见符号，可以通过中英文字符转换(e.g. `:/：, (/（, </《`)后再直接比较。
4. 其他不同种类字符的字符相似度为`0.0`

5. 其他任意字符与`[GAP]`符的相似度均为`0.0`

字符相似度计算的规则是自定义的，目前这种字符相似度计算函数仍需解决的问题是简单忽略了英文字母和汉字读音上的相似性，比如扣和Q、微和V、加和J等。

如果你希望更改或完善字符相似度计算函数，可以直接修改`VariantNW`类中的`__score`函数。



### 2 基于字符相似度匹配的序列对齐算法

**本序列对齐算法以Needleman-Wunsch算法为基础，借鉴其动态规划的思想并引入字符相似度计算函数完成。具体实现细节详见https://jackfromeast.site/2022-02/sequence-alignment-algorithm-based-on-similarity-of-characters.html**

#### 2.1 矩阵初始化

![矩阵初始化](sequence-alignment-algorithm-based-on-similarity-of-characters/矩阵初始化.png)

#### 2.2 矩阵传播

![矩阵传播后](sequence-alignment-algorithm-based-on-similarity-of-characters/矩阵传播后.png)

#### 2.3 回溯

![矩阵回溯](sequence-alignment-algorithm-based-on-similarity-of-characters/矩阵回溯.png)



### 3 运行

```
    variantNW = VariantNW()
   
    # 待对齐的两序列 
    seq1 = "最新辐莉全自动合作睬漂苞蓜日收200--600元 这里→：       hsk65161.com
    seq2 = "最新福利全自动合作彩票包赔日收200--600元这里:hsk65161.com"

    variantNW.set_seqs(seq1, seq2)
    variantNW.propagate()
    
    # 打印NW矩阵
    print(variantNW.NWMatrix)
    
    # 得到对齐后的序列
    aligned_seq1, aligned_seq2 = variantNW.traceback()
    print(aligned_seq1)
    print(aligned_seq2)
    
    # 序列对各位置相似度评分
    print(variantNW.get_aligned_seq_score(aligned_seq1, aligned_seq2))
```
<br>
结果点击放大:
![运行示例](sequence-alignment-algorithm-based-on-similarity-of-characters/%E8%BF%90%E8%A1%8C%E7%A4%BA%E4%BE%8B.png)
