#!/usr/bin/env python
# coding:utf-8

WILDCARD = "."
LEFTKAKKO = "["
RIGHTKAKKO = "]"

from itertools import combinations
import debug
from copy import copy

class IJ(object):
	def __init__(self, i = 0, j = 0):
		self.i = i
		self.j = j

	def __str__(self):
		return "(%s, %s)" %(self.i, self.j)

class Offset_list(object):
	def __init__(self, ijs = []):
		self.ijs = ijs

	def append(self, ij):
		self.ijs.append(ij)

	def cover_strings(self):
		return set([ij.i for ij in self.ijs])

	def __len__(self):
		return len(self.ijs)

	def __str__(self):
		if not self.ijs:
			return "[]"
		return "[" + reduce(lambda x, y:str(x)+str(y), self.ijs) + "]"

class Pattern(object):
	def __init__(self, motif = [], Ls = Offset_list()):
		self.motif = motif
		self.Ls = Ls

	def __str__(self):
		return "".join(self.motif)

	def __eq__(self, other):
		if not isinstance(other, Pattern):
			return False
		if not is_same_list(self.motif, other.motif):
			return False
		return True

	def __len__(self):
		return len(self.motif)

class Elementary_patterns(object):
	"""docstring for Elementary_patterns"""
	def __init__(self, elementary_list = []):
		self.elementary_list = elementary_list

	def insert(self, pattern):
		self.elementary_list.append(pattern)
		# if self.elementary_list:
		# 	for index, elementary_pattern in enumerate(self.elementary_list):
		# 		if not alphabetical(elementary_pattern.motif, pattern.motif):
		# 			self.elementary_list.insert(index, pattern)
		# 			break
		# 	else:
		# 		self.elementary_list.append(pattern)
		# else:
		# 	self.elementary_list.append(pattern)

	def __iter__(self):
		return iter(self.elementary_list)

	def __len__(self):
		return len(self.elementary_list)

	def __getitem__(self, index):
		return self.elementary_list[index]

	def supported(self, config):
		self.elementary_list = [x for x in self.elementary_list if len(x.Ls.ijs) >= config.k]

# @debug.print_in_out
def alphabetical(motif1, motif2):
	# print motif1, "vs", motif2
	length = len(motif2) if len(motif1) > len(motif2) else len(motif1)
	# print length
	for i in range(length):
		if motif1[i] == WILDCARD and motif2[i] == LEFTKAKKO:
			# print "---1"
			# return -1
			return 1
		elif motif1[i] == LEFTKAKKO and motif2[i] == WILDCARD:
			# print "----2"
			# return 1
			return -1
		elif motif1[i] != WILDCARD and motif2[i] == WILDCARD:
			# print "---3"
			# return -1
			return 1
		elif motif1[i] == WILDCARD and motif2[i] != WILDCARD:
			# print "---4"
			# return -1
			return 1
		elif motif1[i] != LEFTKAKKO and motif2[i] == LEFTKAKKO:
			# print "---5"
			# return 1
			return -1
		elif motif1[i] == LEFTKAKKO and motif2[i] != LEFTKAKKO:
			# print "---6"
			# return -1
			return 1
	if len(motif1) >= length:
		# return 1
		return -1
	elif len(motif1) < length:
		# return -1
		return 1
	elif (motif1 < motif2):
		# return 1
		return -1
	else:
		# return -1
		return 1

def is_same_list(list1, list2):
	if len(list1) == len(list2):
		for l1, l2 in zip(list1, list2):
			if not l1 == l2:
				return False
		else:
			return True
	else:
		return False
		

def make_and_insert(substs, eps, current_bit_mask, seqs, i, j):
	temp_pattern = []
	pat_to_add = []
	flag = False
	brack_alter, pat_num, a, b = 0, 0, 0, 0
	temp = ""
	temps = ""

	pat_to_add.append(temp_pattern)
	for w in range(len(current_bit_mask)):
		temp = seqs[i].string[j + w]
		if current_bit_mask[w] == "1":
			for a in range(len(pat_to_add)):
				while len(pat_to_add) < a:
					pat_to_add.append([])
				pat_to_add[a].append(temp)
		elif current_bit_mask[w] == "0":
			for a in range(len(pat_to_add)):
				while len(pat_to_add) < a:
					pat_to_add.append([])
				pat_to_add[a].append(WILDCARD)
		elif current_bit_mask[w] == "2":
			# braなんとかのやつ。あとで実装。
			pass
	ijs = []
	for a in range(len(pat_to_add)):
		count = 0
		# print a, len(eps)
		for it in eps:
			if is_same_list(it.motif, pat_to_add[a]):
				break
			else:
				count += 1
		if count != len(eps):
			# print count, len(eps)
			ij = IJ(i, j)
			# print "---------"
			eps[count].Ls.ijs.append(ij)
		else:
			ij = IJ(i, j)
			ijs.append(ij)
			Ls = Offset_list(ijs)
			eps.insert(Pattern(pat_to_add[a], Ls))
			# print "inserted"
	# print "================================="
	return eps

# def initialize_bit_mask(w, l):
# 	j = w - l
# 	bits = ""
# 	for i in range(l):
# 		bits += "1"
# 	for i in range(j):
# 		bits += "0"
# 	return bits

def perm_bits(n, k):
    result = []
    for bits in combinations(range(n), k):
        s = ['0'] * n
        for bit in bits:
            s[bit] = '1'
        result.append(''.join(s))
    return result

# 1で終わるインデックスを返す
# in	: 001001110100
def search_end_1(bit):
	for i, p in enumerate(bit[::-1]):
		if p == "1":
			return len(bit) - i

def trim0s(bit_mask):
	return bit_mask[:search_end_1(bit_mask)]

def main():
	print alphabetical("W.KK", "M.NA")

if __name__ == "__main__":
	main()
