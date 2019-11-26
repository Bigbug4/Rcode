DIR=getwd()
setwd("E:\\Rcode\\data")

library(dplyr)
data(iris)
attach(iris)
f1 <- filter(iris,Sepal.Length>5,Petal.Length<2)
f1
f2 <- filter(iris,Sepal.Length>5,Petal.Length>2)
f2
f3 <- filter(iris,Sepal.Length<5,Petal.Length<2)
f3
f4 <- filter(iris,Sepal.Length<5,Petal.Length>2)
f4

s1 <- select(iris,Sepal.Length:Petal.Length)[1:5,]
s1
summary(s1)
s1.g <- group_by(s1,Sepal.Width)
summarise(s1,n=n(),l1.m=mean(Sepal.Length),l2.m=mean(Petal.Length))
summarise(s1.g,n=n(),l1.m=mean(Sepal.Length),l2.m=mean(Petal.Length))

a1 <- arrange(iris,Sepal.Length,desc(Petal.Length))[1:5,]
a1
a2 <- arrange(iris,Sepal.Length,desc(Petal.Length)) %>%
       select(Sepal.Length:Petal.Length) %>%
        head(5)
a2

library(ggplot2)
Species.tab <- table(Species)
Species.tab

barplot(Species.tab,ylab = "counts")
qplot(Species,fill=I("#FF0000"),col=I("#000000"),ylab = "counts")

iris$l <- ifelse(Sepal.Length > median(Sepal.Length),"high","low")
iris$l <- as.factor(iris$l)
qplot(Species,data=iris,fill=l,ylab = "counts")
hist(Sepal.Width,breaks = 10,xlab = "width",ylab = "counts")
qplot(Sepal.Width,fill=I(rainbow(13)),col=I("#000000"),ylab = "counts",binwidth=0.2)

boxplot(Petal.Length~Species,col=rainbow(3))
boxplot(Petal.Width~Species,col=rainbow(3))
qplot(Species,Petal.Width,fill=Species,col=I("#000000"),geom="boxplot")
qplot(Species,Petal.Width,fill=I(rainbow(3)),col=I("#000000"),geom="boxplot")

plot(Sepal.Width,Petal.Length,ylab = "Petal.Length")
title(main = "iris",xlab = "",ylab = "")
qplot(Sepal.Width,Petal.Length,col=Species,size=iris$l,
      ylab = "Petal.Length")
qplot(Sepal.Width,Petal.Length,col=iris$l,shape=Species,
      ylab = "Petal.Length")

qplot(Sepal.Width,Petal.Length,col=iris$l,shape=Species,data=iris,
      ylab = "Petal.Length",facets = ~Species)
qplot(Sepal.Width,Petal.Length,col=iris$l,shape=Species,data=iris,
      ylab = "Petal.Length",facets = Species~.)

qplot(Sepal.Width,Petal.Length,col=iris$l,shape=Species,data=iris,
      ylab = "Petal.Length",facets = l~Species)


ggplot(iris,aes(x=Sepal.Width,y=Petal.Length)) + 
  geom_point(aes(shape=Species),color="#FF0000",size=2)
ggplot(iris,aes(x=Sepal.Width,y=Petal.Length)) + 
  geom_point(aes(shape=Species,color=Petal.Width),size=2) 

ggplot(iris) + geom_point(aes(x=Sepal.Width,y=Petal.Length,color=Petal.Width),size=2) 
ggplot(iris) + geom_point(aes(x=Sepal.Width,y=Petal.Length,color=Petal.Width),shape=Species,size=2)
ggplot(iris,aes(x=Sepal.Width,y=Petal.Length)) + 
  geom_line(aes(color=Species),size=1)

ggplot(iris,aes(x=Sepal.Width,y=Petal.Length)) + 
  geom_line(aes(color=Species),size=1) + 
  scale_color_manual(values = c("red", "blue", "green"))

p <- ggplot(iris) + geom_bar(aes(Species,col=I("#000000")),width=0.5,fill="#FF0000") +
  labs(x="Species",title="iris") + theme(plot.title = element_text(hjust = 0.5), panel.grid = element_blank())

p
iris1 <- data.frame(Species.tab) %>% mutate(f=round(Freq/sum(Freq),digits = 3))
ggplot(iris1) + geom_bar(aes(Species,f,col=I("#000000")),width=0.5,fill="#FF0000",stat="identity") +
  scale_y_continuous(labels = scales::percent, limits = c(0,0.35)) +
  xlab("Species") + ylab("Freq") + ggtitle("iris") + theme_bw()

p1<- p + geom_text(aes(Species),stat="count",
            label=paste(iris1$f*100,'%',sep=''),
            colour = "black", vjust=-0.3, size=3)
p1
ggsave("iris2.tiff",dpi = 500,plot = p1)

ggplot(iris) + geom_histogram(aes(Petal.Length,col=I("#000000")),fill="#FF0000",binwidth=0.2)
ggplot(iris) + geom_density(aes(Petal.Length,col=I("#000000")),fill="#FF0000")
ggplot(iris,aes(x=Sepal.Width,y=Petal.Length)) + 
  geom_area(aes(fill=Species),col="#000000")

ggplot(iris,aes(x=Sepal.Width,y=Petal.Length)) + 
  geom_point(aes(color=Species,shape=Species),size=2) + 
  labs(title = "Differential Species",x="Sepal.Width",y="Petal.Length") +
  theme(legend.position="right", legend.title = element_blank(),panel.grid = element_blank())

ggsave("iris1.tiff",dpi = 500)


