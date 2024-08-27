
# Due to Giacomo Pope
def random_point(E):
	"""
	Returns a random point on the elliptic curve E
	assumed to be in Montgomery form with a base 
	field which characteristic p = 3 mod 4
	"""
	A = E.a_invariants()[1]
	C = E.a_invariants()[3]
	if E.a_invariants() != (0,A,0,C,0):
		raise ValueError("Function `random_point` assumes the curve E is in the Montgomery model")

	# Try 10000 times then give up, just protection
	# for infinite loops
	F = E.base_ring()
	for _ in range(10000):
		x = F.random_element()
		y2 = x*(x**2 + A*x + C)
		if y2.is_square():
			y=y2.sqrt()
			return E(x,y)

	raise ValueError("Generated 10000 points, something is probably going wrong somewhere.")

def int_to_digit_t(x,n_words):
	L_words=[]
	y=x
	for i in range(n_words):
		L_words.append(hex(y%(2**64)))
		y=y//(2**64)
	return L_words
