# Analyzed bootstrap Exit time & survival

rm(list = ls())
graphics.off()

# Analyze bootstrapped exit times

load(file='ET+Sboot_2021.Rdata') #--------------------------------------------------------------
LET03 = ETLR[,1]  # convert from 5 min to hours
RET03 = ETLR[,2]

# Nominal exit times from the data
SmeanETl = 7.44 # hours
SmeanETr = 13.76 # hours

# Correct bootstrap bias
medLET = median(LET03)
medRET = median(RET03)
#
Lbias = SmeanETl - medLET
Rbias = SmeanETr - medRET
#
LET30 = LET03 + Lbias
LET3 = subset(LET30,subset=(LET30 > 1)) # remove negative ET due to bias correction
RET30 = RET03 + Rbias
RET3 = subset(RET30,subset=(RET30 > 1))

probseq = c(0.05,0.1,0.25,0.5,0.75,0.9,0.95)
print('',quote=F)
print('2021 Exit Time ============================================',quote=F)
print('Left basin',quote=F)
print('quantiles',quote=F)
print(quantile(LET3,probs=probseq),quote=F)
print('summary',quote=F)
print(summary(LET3),quote=F)
print(c('sd = ',sd(LET3)),quote=F)

print('---------------------------------',quote=F)
print('Right basin',quote=F)
print('quantiles',quote=F)
print(quantile(RET3,probs=probseq),quote=F)
print('summary',quote=F)
print(summary(RET3),quote=F)
print(c('sd = ',sd(RET3)),quote=F)

# Analyze bootstrapped half life =========================================
#
load(file='ET+Sboot_2021.Rdata') #--------------------------------------------------------------
LET03 = ShalfLR[,1]  # convert from 5 min to hours
RET03 = ShalfLR[,2]

# Nominal half life from the data
SmeanL = 4.28 # hours
SmeanR = 40.0 # hours

# Correct bootstrap bias
medLET = median(LET03)
medRET = median(RET03)
#
Lbias = SmeanL - medLET
Rbias = SmeanR - medRET
#
LET30 = LET03 + Lbias
LET3 = subset(LET30,subset=(LET30 > 1)) 
RET3 = RET03 + Rbias

probseq = c(0.05,0.1,0.25,0.5,0.75,0.9,0.95)

print('',quote=F)
print('2021 Half Life =========================================',quote=F)
print('Left basin',quote=F)
print('quantiles',quote=F)
print(quantile(LET3,probs=probseq),quote=F)
print('summary',quote=F)
print(summary(LET3),quote=F)
print(c('sd = ',sd(LET3)),quote=F)

print('---------------------------------',quote=F)
print('Right basin',quote=F)
print('quantiles',quote=F)
print(quantile(RET3,probs=probseq),quote=F)
print('summary',quote=F)
print(summary(RET3),quote=F)
print(c('sd = ',sd(RET3)),quote=F)

