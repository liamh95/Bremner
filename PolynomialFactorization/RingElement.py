class RingElement:

	def __init__(self, ring):
		self.ring = ring

	def __add__(self, other):
		assert self.ring == other.ring
		