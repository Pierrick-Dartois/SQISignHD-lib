from sage.all import *

class Tree:
	def __init__(self,node):
		self._node=node
		self._edges=[]
		self._children=[]

	def add_child(self,child,edge):
		self._children.append(child)
		self._edges.append(edge)

	def look_node(self,node):
		if self._node==node:
			return self
		elif len(self._children)>0:
			for child in self._children:
				t_node=child.look_node(node)
				if t_node!=None:
					return t_node

	def edge_product(self,L_factors,factor_node=ZZ(1)):
		n=len(self._children)
		L_prod=[(factor_node,self._node)]
		for i in range(n):
			L_prod+=self._children[i].edge_product(L_factors,factor_node*L_factors[self._edges[i]])
		return L_prod

