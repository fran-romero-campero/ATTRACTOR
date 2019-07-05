ggplot(network.data, aes(x.pos,y.pos)) + 
  theme(panel.background = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks.y = element_blank()) + 
  geom_point(color=node.colors,size=14)


mean.expression[1:4,1:2]

min(mean.expression,na.rm = T)
max(mean.expression,na.rm=T)

help("spline")
plot(spline(x = seq(from=0,to=24,by=4),
       y = 3 + as.vector(scale(c(mean.expression[7,],mean.expression[7,1]))),n = 48),type="o")


normalized.interpolated.expression.data <- matrix(0,nrow= nrow(mean.expression), ncol=48 )
  
expression.data <- mean.expression
expression.data[is.na(expression.data)] <- 1
expression.data[1,]
i <- 1
expression.data[i,]
for(i in 1:nrow(mean.expression))
{
  normalized.interpolated.expression.data[i,] <- spline(x = seq(from=0,to=24,by=4),
         y = 3 + as.vector(scale(c(expression.data[i,],expression.data[i,1]))),n = 48)$y
}



head(normalized.interpolated.expression.data)




help(scale)

plot((as.vector(scale(mean.expression[5,]))),type="o")

plot(mean.expression[5,]/(max(mean.expression[5,])),type="o")





