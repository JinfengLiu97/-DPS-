# 在JJ过程中测量微分截面、其系统误差、绘图并测量DPS占比  

四维拟合需要较高版本的CMSSW，否则会报错。目前确定的是CMSSW_12_4_0是可以工作的。在四维拟合之前要运行cmsenv。  

这里只介绍代码的运行，其底层原理需要自己去看AN。  

## 样本  

所有的样本都是过筛后添加了接受度&效率权重（Weight_sum）的样本。

`Data_MixWeighted.root`: 使用SPS/DPS(0.3/0.7)混合接受度&效率进行修正的样本。提供主要结果；  
`Data_SPSWeighted.root`: 使用SPS接受度&效率进行修正的样本。用来估计修正带来的误差；  
`Data_DPSWeighted.root`: 使用DPS接受度&效率进行修正的样本。用来估计修正带来的误差；  

`MC_SPS.root`: SPS的MC样本；  
`MC_DPS.root`: DPS的MC样本。  

## 计算微分截面  

此处以双J不变质量（`JJMass`，变量名为`FourMuonMass_`）为例。  

计算微分截面就是在不同的区间（Bin）中重新进行拟合获得信号产额，所以核心代码还是四维拟合的代码（`Fit_Final.C`）。但是为了便于连续计算，所以添加了一些接口。现在主函数的形式为：  
`void Fit_Final(float FourMuonMass_LowCut_Value, float FourMuonMass_HighCut_Value, int MassNumber, TString Output_Name) `  
参数分别是：目标区间的下限；目标区间的上限；目标区间的编号（一般由小至大从1开始）；输出文件的文件名  

使用脚本`AutoRun.C`来连续调用主函数。脚本中设置了每个区间的上下限。当需要调整分隔时可以在这里修改。`AutoRun.C`中还规定了输出文件的文件名（Record）。  
在`Differential/JJMass`中运行：  
`root -l AutoRun.C > Log.log`  
来计算微分截面（会产生大量日志信息，建议重定向，否则屏幕会被淹没）。运行成功的话会报告生成了一系列图片文件，并产生一个输出文件`Record.txt`  

需要检查`Record.txt`。其中会显示每个区间的：编号，[下限，上限]， 拟合状态(status)，以及各个成分的产额及其统计误差。需要确保每个区间的拟合都是成功的(status=0)，且产额及误差是合理的（如果产额都是拟合时设定的初值，误差极小就是不正常）。还可以检查Plot目录下的四维拟合图片来确定拟合情况。  

如果出现部分区间拟合状况不佳的情况，可以单独调用Fit_Final函数：  
运行：  
`root -l Fit_Final.C > Log.log`
再重新调用主函数，比如：  
`Fit_Final(7.5, 17.5, 1, "Record")`  
来重新对目标区间进行拟合（与此同时需要双手合十真诚祈祷）。  

需要注意的是，如果不进行额外的设置（图片输出的控制在`Fit_Final.C`的后半段），那么新产生的图片会直接覆盖之前的图片。但是输出文件（`Record.txt`）的写入模式是附加模式（`a`）。所以如果前后的输出文件是同名的，那么新写入的内容会追加到文件的末尾，不会产生覆盖。  

可以修改`AutoRun.C`与`Fit_Final.C`中的变量名等参数来调整测量的维度。可以参考数据样本中提供的branch。本研究中主要测量的维度及其区间是：  
双J的不变质量：`FourMuonMass_` [7.5, 17.5, 27.5, 37.5, 47.5, 57.5, 67.5, 77.5]以及一个overflow区间；  
两个J的快度差：`DeltaY_` [0.0, 0.5, 1.0, 1.5, 2.0, 2.5]以及一个overflow区间；  
双J的横动量：`FourMuonPt_` [0, 5, 10, 15, 20, 25, 30, 35, 40]以及一个overflow区间；  
双J的快度：`FourMuonY_` [0, 0.4, 0.8, 1.2, 1.6, 2.0]；  
两个J的方位角差：`DeltaPhi_` [0, pi/8, pi/4, 3pi/8, pi/2, 5pi/8, 3pi/4, 7pi/8, pi]；  
两个J的横动量的平均值：`JPt` [10, 15, 20, 25, 30, 40]。  

## 计算系统误差  

仍然以双J不变质量（`JJMass`，变量名为`FourMuonMass_`）为例。  

系统误差就是在不同的特殊条件下重复测量微分截面，所以核心代码还是四维拟合代码（`Fit_Final.C`）。但是额外提供了一个参数：`Standard_Value`，通过其向函数传递该区间产额的标准值（即上一步得到的信号产额）来直接计算系统误差。相应的，`AutoRun.C`也添加了相应的数组来传递标准值。另外由于要计算不同的系统误差，所以不同的拟合代码也进行了相应的特化处理。  

具体的对系统误差的说明参考AN。  

系统误差需要逐项计算，再手动计算总误差

### 亮度的误差  

所有区间都是0.120%    

### 衰变分支比的误差  

所有区间都是1.106%  

### 修正带来的误差（`Systematic/Correction`）  

需要换用不同修正的样本来重新对微分截面进行测量  

先估计单纯使用SPS修正造成的误差。  

在`Systematic/Correction/JJMass/SPS`中，如果有需要，可以修改`AutoRun.C`中的区间分隔，标准值与输出文件。运行:  
`root -l AutoRun.C > Log.log`  
此处调用的`Fit_Final.C`读入的数据文件是单纯使用SPS接受度&效率进行修正的样本（`Data_SPSWeighted.root`）。与微分截面的计算类似的，会产生一个输出文件（`SPS_Record.txt`）并在Plot目录下产生一系列四维拟合图片。输出文件与计算微分截面时的输出文件类似，但是在每个区间的末尾给出了系统误差的数值（没有单位，也即需要自己把小数点向后挪两位得到百分比的系统误差）。  

计算其他系统误差时产生的输出文件和图片的内容大都与之类似，以下不再赘述。  

仍然需要检查每个区间，确保拟合正常。单独对某个区间进行拟合的方法和微分截面的计算一样，但是要额外给出该区间信号产额的标准值，如：
`Fit_Final(7.5, 17.5, 1, 2858, "SPS_Record")`  
（别忘记祈祷哦）  

类似的，再估计单纯使用DPS修正造成的误差。在`Systematic/Correction/JJMass/DPS`中运行：  
`root -l AutoRun.C > Log.log`  
此处调用的`Fit_Final.C`读入的数据文件是单纯使用DPS接受度&效率进行修正的样本（`Data_DPSWeighted.root`）。  

修正带来的误差需要比较使用SPS修正造成的误差和使用DPS修正造成的误差，取其中较大的那个。  

***自动化可以进一步提高。可以在一个`Fit_Final.C`中读取两个数据文件（`Data_SPSWeighted.root`与`Data_DPSWeighted.root`），进行两次拟合，产生两次输出，并直接比较两次拟合误差的最大值。这样就避免了手动比较两者误差的麻烦。但是当某些区间拟合情况不佳时，处理时会产生一次额外的拟合。另外，画图的代码也需要额外的处理。**  

### ctau的分布带来的系统误差（`Systematic/CtauShape`）  

需要换用不同的瞬发Jpsi的ctau分布参数来重新对微分截面进行测量。  

先估计单纯使用SPS的ctau分布造成的误差。在`Systematic/CtauShape/JJMass/SPS`中运行：  
`root -l AutoRun.C > Log.log`  
此处调用的`Fit_Final.C`使用的瞬发Jpsi的ctau分布的参数来自于之前对一个SPS样本的拟合。  

再估计单纯使用SPS的ctau分布造成的误差。在`Systematic/CtauShape/JJMass/DPS`中运行：  
`root -l AutoRun.C > Log.log`  
此处调用的`Fit_Final.C`使用的瞬发Jpsi的ctau分布的参数来自于之前对一个DPS样本的拟合。  

ctau的分布带来的系统误差需要比较使用SPS的ctau分布造成的误差和使用DPS的ctau分布造成的误差，取其中较大的那个。  

***自动化的提高与修正的类似。可以在`Fit_Final.C`中构造两个ctau分布，再进一步构造两个拟合的PDF，进行两次拟合。**  

### 拟合子稳定性带来的系统误差（`Systematic/Fitter`）  

在`Systematic/Fitter/JJMass`中运行：  
`root -l AutoRun.C > Log.log`  
此处调用的`Fit_Final.C`会进行多次拟合来构建一个对拟合子稳定性的测试。首先对（目标区间内的）数据样本进行一次拟合，假设拟合得到的产额是：`PP(2000)`, `PNP(1000)`, `NPNP(12000)`, `JpsiMuMu(800)`, `MuMuMuMu(20)`  
随后代码会（重复五次）产生一定量的toy MC事例，附加到数据样本上再次拟合：  
第一次：在数据样本上附加2000个PP事例，再进行拟合（预期得到的各成分产额是`PP(4000)`, `PNP(1000)`, `NPNP(12000)`, `JpsiMuMu(800)`, `MuMuMuMu(20)`）；  
第二次：在数据样本上附加1000个PNP事例，再进行拟合（预期得到的各成分产额是`PP(2000)`, `PNP(2000)`, `NPNP(12000)`, `JpsiMuMu(800)`, `MuMuMuMu(20)`）；  
....  
第五次：在数据样本上附加20个MuMuMuMu事例，再进行拟合（预期得到的各成分产额是`PP(2000)`, `PNP(1000)`, `NPNP(12000)`, `JpsiMuMu(800)`, `MuMuMuMu(40)`）。  
每次拟合，代码都会计算信号产额的误差（实际拟合结果与预期结果之间的误差），并自动给出五次拟合中最大的误差。  
具体的拟合结果可以在输出文件`Record.txt`中查看。其中每个区间都会给出五次拟合的拟合情况、附加的种类、各成分产额与信号产额的相对误差，并在该区间的所有拟合结束后给出该区间不同拟合中的最大误差。以该最大误差作为拟合子稳定性带来的误差。  

另外，由于进行了大量的拟合，所以我直接注释掉了`Fit_Final.C`中绘图的代码，也即不会产生任何图片文件。如果希望查看图片，需要仔细修改绘图代码，确保每个区间内不同拟合的图片不会产生覆盖。  

### 参数固定带来的系统误差（`Systematic/Fix`）  

需要依次松动PDF2/PDF3/PDF4（具体参考AN）的形状参数来重新对微分截面进行测量。  

先估计松动PDF2的参数造成的误差。在`Systematic/Fix/JJMass/PDF2`中运行：  
`root -l AutoRun.C > Log.log`  
此处调用的`Fit_Final.C`松动了PDF2所有的形状参数。  

类似的，依次在`Systematic/Fix/JJMass/PDF3`和`Systematic/Fix/JJMass/PDF4`中运行`AutoRun.C`。  

比较三次拟合中的误差，取最大的那个作为参数固定造成的误差。  

***自动化的提高与ctau形状类似，可以构造三个不同的拟合PDF来进行三次拟合。**  

### 寿命变量选择带来的系统误差（`Systematic/SigLxy`）  

在`Systematic/SigLxy/JJMass`中运行：  
`root -l AutoRun.C > Log.log`  
此处调用的`Fit_Final.C`使用了SigLxy（LxyPV的显著度，参考AN）代替了ctau。  

在计算这一项误差的时候需要注意，由于未知原因，该拟合十分耗时。每次拟合可能需要十数分钟。建议找个网好的地方，打开一堆窗口，每个窗口运行一个维度，然后就可以去睡觉了。  

### 总系统误差  

当完成所有系统误差的计算后，就可以手动计算总系统误差。将各项系统误差平方后求和，再开根，就得到了总系统误差（**别忘记亮度和分支比的系统误差**）。  

## 图片绘制  

在完成所有误差的计算后，就可以开始绘制微分截面的分布（`Painter/DifferentialPainter.C`）和系统误差的分布（`Painter/SystematicPainter.C`）。  

在Painter中修改`DifferentialPainter.C`。我已经填入了一组示例性的数字，你可以修改：  

区间分隔：XXBinEdge （如果修改了，要相应修改后续数组的长度（紧跟在区间分隔代码后的四个数组）、填充时的循环上限（数组之后的一行、后续的条形图填充内容以及输出部分）、声明条形图时的区间数量、计算微分截面时的区间宽度、绘图时横轴显示的范围以及overflow字样的显示位置）；  
各区间信号产额：XXBinContent；  
各区间信号产额统计误差：XXBinError_Statistic （单位为1）；  
各个区间的总系统误差：XXBinError_Systematic （单位为%）。 

运行：  
`root -l DifferentialPainter.C`  
会产生相应的微分截面分布图片（`Painter/DifferentialPlot/`，pdf格式）。可以打开图片并进一步调整代码。
同时，运行该代码会直接输出所有微分截面、其统计误差、其系统误差的具体数值。已经修改了它的格式，使之符合latex的排版，这样可以方便的编辑AN中的表格（但是由于表格中不同截面的有效数字现在取得不太一样，所以可能还需要手动调整有效数字）。  

在Painter中修改`SystematicPainter.C`。类似的，示例性的数字已经填入，你可以修改：  

区间分隔：XXBinEdge（如果修改了，要相应修改后续数组的长度（在构建Graph时的三个数组）、初始化数组时的循环上限、填充时的循环上限（数组后的循环）、声明Graph时的区间数量以及绘图时横轴显示的范围）；  
各项系统误差的数值 （单位为%）。  

运行：  
`root -l SystematicPainter.C`  
会产生相应的系统误差的分布（`Painter/SystematicPlot/`，png格式（我觉得png格式的好看一些（这个图片又不会正式发表）。可以打开图片并进一步调整代码。  

***自动化可以进一步提高，或许可以修整之前的输出文件的格式，使之可以被容易复制到绘图代码中。也可以在绘图代码中添加读入文件直接读入之前产生的输出文件（画系统误差时一个维度就要手动读几十个数字，太折磨了）。如果自动读入可以实现的话，系统误差的绘图或许可以每个维度只构建一个Graph，然后重复对其进行填入并绘图的程序。**  

## fDPS的测量  

测量完微分截面后就可以测量DPS的占比了（`fDPS/`）。主要方法是在DeltaY和DeltaPhi维度上进行同步模板拟合。  

在`fDPS`中修改`SimFit.C`，需要将DeltaY和DeltaPhi的微分截面信息输入。我已经填入了一组示例性数据。数据内容和微分截面的绘图代码中（`Painter/DifferentialPainter.C`）的相应内容是完全一样的。代码会自动读入相应的MC（`Sample/MC_SPS.root与Sample/MC_DPS.root`），并加载拟合使用的模板。  

运行：  
`root -l SimFit.C`
会产生相应的拟合图片（`fDPS/Plot/`）。可以打开图片对拟合进行检查。该代码会在最后输出fDPS的测量值及其统计误差。

***现在这个代码的主要问题是，我暂时没有找到将DeltaY的overflow区间加入到拟合中的方法，需要进一步改进。**
