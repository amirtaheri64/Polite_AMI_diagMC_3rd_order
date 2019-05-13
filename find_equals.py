

a=[9, 16, 26,35,14,33, 19,28, 28,29, 33,34,12,26,9, 32]  # Array of equals

a_pairs=[]
for i in range (0,len(a)/2):
  a_pairs.append([a[2*i],a[2*i+1]])
print a_pairs
c=[]
if len(a)!=0:
  c.append(a[0])
  c.append(a[1])

for j in range (0,len(a_pairs)):
  for i in range(0,len(a_pairs)):
    if (a_pairs[i][0] in c or a_pairs[i][1] in c):
      if a_pairs[i][0] not in c:
        c.append(a_pairs[i][0])
      if a_pairs[i][1] not in c:
        c.append(a_pairs[i][1])
  for k in range (0,len(a)):
    if a[k] not in c:
      c.append(0)
      c.append(a[k])
      break
      
print c  # Array of equal islands
