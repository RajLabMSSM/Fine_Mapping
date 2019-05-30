library(ggplot2)
library(data.table)
x = read.table(
  "/Users/evan/Desktop/MyND/eQTL/boxplot.tmp.txt",
  header = T,
  row.names = 1,
  check.names = F
)

y = data.frame(t(x), check.names = F)
head(y)
colnames(y) = gsub(':', '', colnames(y))
y$rs814228 = as.factor(y$chr1240383902)
y$rs35355673 = as.factor(y$chr1488475921)
y$rs1964536 = as.factor(y$chr723380185)


p = ggplot(y, aes(rs1964536, GPNMB, color = rs1964536))
p + geom_boxplot(width = 0.4) + geom_jitter(position = position_jitter(w = 0.1, h = 0), alpha = 0.9) +
  theme_bw() + ylab("Expression (Residuals)") + xlab("chr7:23380185 \n pGWAS = 2.30E-13, MAF = 0.447") +
  labs(title = "rs1964536 (GPNMB)", subtitle = "P-Value = 4.31E-02, b = -0.285") +
  scale_x_discrete(labels = c("0" = "TT", "1" = "TC",
                              "2" = "CC")) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
    plot.subtitle = element_text(
      size = 18,
      hjust = 0.5,
      color = "black"
    ),
    text = element_text(size = 22),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(size = 22)
  ) + scale_fill_brewer(palette = "Accent")
