from sage.all import *
from utilities.strategy import precompute_strategy_with_first_eval, precompute_strategy_with_first_eval_and_splitting
from isogenies.isogeny_dim4 import IsogenyDim4


class IsogenyChainDim4:
	def __init__(self, B_K, first_isogenies, e, m, splitting=True, strategy = None):
		self.e=e
		self.m=m

		if strategy == None:
			strategy = self.get_strategy(splitting)
		self.strategy = strategy

		self._isogenies=self.isogeny_chain(B_K, first_isogenies)


	def get_strategy(self,splitting):
		if splitting:
			return precompute_strategy_with_first_eval_and_splitting(self.e,self.m,M=1,S=0.8,I=100)
		else:
			return precompute_strategy_with_first_eval(self.e,self.m,M=1,S=0.8,I=100)

	def isogeny_chain(self, B_K, first_isogenies):
		"""
		Compute the isogeny chain and store intermediate isogenies for evaluation
		"""
		# Store chain of (2,2)-isogenies
		isogeny_chain = []

		# Bookkeeping for optimal strategy
		strat_idx = 0
		level = [0]
		ker = B_K
		kernel_elements = [ker]

		# Length of the chain
		n=self.e-self.m
		
		for k in range(n):
			prev = sum(level)
			ker = kernel_elements[-1]

			while prev != (n - 1 - k):
				level.append(self.strategy[strat_idx])
				prev += self.strategy[strat_idx]

				# Perform the doublings and update kernel elements
				# Prevent the last unnecessary doublings for first isogeny computation
				if k>0 or prev!=n-1:
					ker = [ker[i].double_iter(self.strategy[strat_idx]) for i in range(4)]
					kernel_elements.append(ker)

				# Update bookkeeping variable
				strat_idx += 1

			# Compute the codomain from the 8-torsion
			if k==0:
				phi = first_isogenies
			else:
				phi = IsogenyDim4(Th,ker)

			# Update the chain of isogenies
			Th = phi._codomain
			isogeny_chain.append(phi)

			# Remove elements from list
			if k>0:
				kernel_elements.pop()
			level.pop()

			# Push through points for the next step
			kernel_elements = [[phi(T) for T in kernel] for kernel in kernel_elements]

		return isogeny_chain

	def evaluate_isogeny(self,P):
		Q=P
		for f in self._isogenies:
			Q=f(Q)
		return Q

	def __call__(self,P):
		return self.evaluate_isogeny(P)

	def dual(self):
		n=len(self._isogenies)
		isogenies=[]
		for i in range(n):
			isogenies.append(self._isogenies[n-1-i].dual())
		return DualIsogenyChainDim4(isogenies)


class DualIsogenyChainDim4:
	def __init__(self,isogenies):
		self._isogenies=isogenies

	def evaluate_isogeny(self,P):
		n=len(self._isogenies)
		Q=P
		for j in range(n):
			Q=self._isogenies[j](Q)
		return Q

	def __call__(self,P):
		return self.evaluate_isogeny(P)




		

