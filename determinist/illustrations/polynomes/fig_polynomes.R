library("ggplot2")
library("pracma")


s = 0.25; c = 0.25 
#s = 0.35; c = 0.25
#s = 0.75; c = 0.25
df_h <- data.frame(h = linspace(0.01,0.99,99),
                   p_star_zyg = (s*(1-(1-c)*(1-linspace(0.01,0.99,99))) - c*(1-s))/(s*(1-2*(1-c)*(1-linspace(0.01,0.99,99)))),
                   p_star_ger = ((1-s*linspace(0.01,0.99,99))*(1+c)-1)/(s*(1-2*linspace(0.01,0.99,99))),
                   ch_line = (1-c)*(1-linspace(0.01,0.99,99)))

ggplot(data=df_h, aes(x=h)) + 
  geom_line(aes(y = p_star_zyg), color = "darkred") + 
  geom_line(aes(y = p_star_ger), color="steelblue") +
  geom_line(aes(y = ch_line), color="grey") +
  geom_vline(xintercept=0, size=0.1) + geom_hline(yintercept=0, size=0.1) +
  scale_y_continuous(name="f(p)", limits=c(-1.1, 2.1)) +
  theme(axis.text = element_text(size = 5),text = element_text(size = 8, family = "Times New Roman")) 





s = 0.75  # when selection only acts on survival with a fitness cost s 
h = 0  # and sh for heterozygous individuals
c = 0.5              # homing rate
homing = "zygote"   # "zygote" or "germline"
num = "8"
label_size = 17

s_1 = c/(1-h*(1-c))  
if(homing == "zygote")
  {s_2 = c/(2*c + h*(1-c))
  A_str = s*(2*(1-c)*(1-h)-1)
  B_str = s*((1-(1-c)*(1-h)))
  if ((1-c)*(1-h)!= 1/2){p_star = (s*(1-(1-c)*(1-h)) - c*(1-s))/(s*(1-2*(1-c)*(1-h)))}
  if(s_1 < s_2)
    {cat("\nCoexistence", ":", round(s_1,3), "< s <", round(s_2,3))
  }else{
    cat("\nBistability", ":", round(s_2,3), "< s <", round(s_1,3))}
  if(c*(1-2*s)-(1-c)*s*h > 0){cat("Linear speed :", 2*sqrt(c*(1-2*s)-(1-c)*s*h))} 
  }

if(homing == "germline")         
  {s_2 = c/(2*c*h + h*(1-c))
  A_str = s*(1-2*h)
  B_str = s*h
  if (h!= 1/2){p_star = ((1-s*h)*(1+c)-1)/(s*(1-2*h))}
  if(s_1 < s_2){cat("\nCoexistence", ":", round(s_1,3) , "< s <", round(s_2,3))
  }else{cat("\nBistability", ":", round(s_2,3), "< s <", round(s_1,3))}
  if(c*(1-2*s*h)-(1-c)*s*h > 0){cat("Linear speed :", 2*sqrt(c*(1-2*s*h)-(1-c)*s*h))}
  }


df <- data.frame(p=linspace(-1,2,300),
                 f=((-A_str*linspace(-1,2,300) - B_str + c*(1-s))*(1-linspace(-1,2,300))*linspace(-1,2,300))/(-A_str*linspace(-1,2,300)**2 - 2*B_str*linspace(-1,2,300) + 1))
if((homing=="zygote"& (1-c)*(1-h)==1/2) | (homing=="germline"& h==1/2)){df_points <- data.frame(x=c(0,1),y=c(0,0), col=c("black", "black"))
}else{df_points <- data.frame(x=c(0,1,p_star),y=c(0,0,0), col=c("black", "black", "red"))}

# Graphique linéaire simple avec des points
ggplot(data=df, aes(x=p, y=f, group=1)) +
  geom_line(color="steelblue") + 
  geom_vline(xintercept=0, size=0.1) + geom_hline(yintercept=0, size=0.1) + 
  geom_point(data=df_points, aes(x=x, y=y, group=1, colour = col), shape=18, size=2) + scale_colour_identity() +
  scale_x_continuous(name="p", limits=c(-1, 2)) + 
  scale_y_continuous(name="f(p)", limits=c(-0.2, 0.2)) +
  theme(axis.text = element_text(size = 5),text = element_text(size = 8, family = "Times New Roman")) 
ggsave(paste0("Desktop/article/courbes_r/", homing, num,".svg"), plot = last_plot(), device = "svg")




