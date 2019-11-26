
tiff(filename = "ex1.tiff",pointsize = 20)
par(mfrow=c(1,3))

plot(sin, -pi, 2*pi,type = 'l',col="#FF0000",xlab='x',ylab='y')
plot(cos, -pi, 2*pi,type = 'l',col="#00FF00",xlab='x',ylab='y')
plot(sin, -pi, pi,type = 'l',col="#0000FF",xlab='x',ylab='y')
dev.off()