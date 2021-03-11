library(basicTrendline)
x = c(3.6,4.6,4.7,5.7,6.6,7.5,8.3,9.4,10.1,11.12,13.9,13.5,14.2,13.7,15.6)
y = c(11,11,12,13,14,15,16,17,18,19,20,21,22,23,25)
#自动添加95%置信区间lines and fill color
trendline(x, y, model="line2P", ePos.x = "topleft", summary=TRUE, eDigit=5)
#只添加95%置信区间的lines，不fill color 
#(set CI.fill = FALSE)
trendline(x, y, model="line3P", CI.fill = FALSE, CI.color = "black", CI.lty = 2, linecolor = "blue")
#只绘制回归曲线，不添加95%置信区间 
#(set CI.color = NA)
trendline(x, y, model="log2P", ePos.x= "top", linecolor = "red", CI.color = NA)
#显示方程，不显示R值和P值 
#(set show.Rpvalue = FALSE)
trendline(x, y, model="line3P", show.equation = TRUE, show.Rpvalue = FALSE)
#自定义方程中的参数的名称:‘xname’, ‘yname’, ‘yhat’, ‘Rname’, ‘Pname’
trendline(x, y, model="line3P", xname="a", yname=paste(beta^15,b),yhat=FALSE, Rname=1, Pname=0, ePos.x = "bottom")
#改变方程的 小数位，字体颜色，字号大小
trendline(x, y, model="power2P", ePos.x = "topleft", summary=TRUE,eDigit = 3, eSize = 1.4, text.col = "blue")
#不显示方程，只显示回归曲线 
#(set ePos.x = NA)trendline(x, y, model="power2P",ePos.x = NA)
#设置绘图区大小
par(mgp=c(1.5,0.4,0), mar=c(3,3,1,1), tck=-0.01, cex.axis=0.9)
trendline(x, y)
