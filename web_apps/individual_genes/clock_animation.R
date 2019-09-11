#Function for radian conversion
radian.conversion <- function(alpha)
{
  rad <- (alpha*pi/180)
  return(rad)
}

## Generate coordinates for inner and outer circle in the clock representation for 
## clock visualizer
angle <- seq(from=0, to=2*pi, by=0.01)
x.circle.1 <- radius.1*sin(angle)
y.circle.1 <- radius.1*cos(angle)

radius.2 <- radius.1 - radius.1/12
x.circle.2 <- radius.2 * sin(angle)
y.circle.2 <- radius.2 * cos(angle)


#Plot circle
par(mar=c(0,0,0,0))
plot(x.circle.1,y.circle.1, type = "l", lwd=3, axes=FALSE, xlab = "", ylab="",xlim=c(-1.2 * radius.1, 1.2 * radius.1),ylim=c(-1.2 * radius.1, 1.2 * radius.1))
lines(x.circle.2, y.circle.2, lwd=3)
x.polygon <- c(sin(seq(from=0, to=-pi, by=-0.01)) * radius.2, 
               sin(seq(from=-pi, to=0, by=0.01))* radius.1)
y.polygon <-c(cos(seq(from=0, to=-pi, by=-0.01)) * radius.2, 
              cos(seq(from=-pi, to=0, by=0.01))*radius.1)
polygon(x = x.polygon, y = y.polygon, col = "black")
for (j in 0:5)
{
  angle.zt <- radian.conversion(alpha = 60*j)
  zt <- 4*j
  current.zt <- paste("ZT", zt,  sep = "")
  text(x = (radius.1 + radius.1/6)*sin(angle.zt), y = (radius.1 + radius.1/6)*cos(angle.zt), labels = current.zt,cex = 1.5,font=2)
  lines(x = c(radius.1 * sin(angle.zt), (radius.1 + radius.1/20)* sin(angle.zt)), 
        y = c(radius.1 * cos(angle.zt), (radius.1 + radius.1/20)* cos(angle.zt)), lwd=2)
}

radio.flecha <- 80

i <- 3
angle.zt <- radian.conversion(alpha = 60*i)
angle.zt.sec <- radian.conversion(alpha = 6*i) ##10 veces más rápido, 10 veces más pequeño el ángulo

# lines(x = c(0, sin(angle.zt)*radio.flecha), y = c(0, cos(angle.zt)*radio.flecha), lwd = 5)
arrows(x0 = 0, y0 = 0, x1 = sin(angle.zt)*radio.flecha, y1 = cos(angle.zt)*radio.flecha,lwd = 5)
arrows(x0 = 0, y0 = 0, x1 = sin(angle.zt.sec)*radio.flecha, y1 = cos(angle.zt.sec)*radio.flecha,lwd = 5)

