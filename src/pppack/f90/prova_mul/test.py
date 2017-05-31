import copy as cp
class P(object):
  def __init__(self,k):
    try:
      self.k=[j for j in k]
    except TypeError:
      self.k=[k]
  def __mul__(self,P0):
    return P(cp.copy(self.k)+cp.copy(P0.k))
  def __str__(self): return str(self.k)

P1,P2=P(1),P(3)
P3=P1*P2

print P1,P2,P3

P1.k=5
print P1,P3
