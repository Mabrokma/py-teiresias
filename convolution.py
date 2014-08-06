#!/usr/bin/env python
# coding:utf-8

from pattern import *


def remove_entries(p, dir_p, dir_s, overlap, delim):
	dir_p.delete(p, overlap, delim)
	dir_s.delete(p, overlap, delim)

##############################################################
##					Convolute modules						##
##############################################################

def right_convolute(t, pQ, suffix_w):
	pat = Pattern([], Offset_list([]))
	new_pattern = t.motif + pQ.motif[len(suffix_w):]
	pat.motif = new_pattern
	delta_j = non_bracket_diff_length(t.motif, suffix_w)
	pol = match_ijs_rhs(t, pQ, delta_j)
	pat.Ls = pol
	return pat

def match_ijs_rhs(t, pQ, delta_j):
	pol = Offset_list([])
	for itPIJ in t.Ls.ijs:
		for itSIJ in pQ.Ls.ijs:
			match_ijs(itPIJ, itSIJ, delta_j, pol)
	return pol

def left_convolute(pQ, t, prefix_w):
	pat = Pattern([], Offset_list([]))
	new_pattern = pQ.motif + t.motif[len(prefix_w):]
	pat.motif = new_pattern
	delta_j = non_bracket_diff_length(pQ.motif, prefix_w)
	pol = match_ijs_lhs(pQ, t, delta_j)
	pat.Ls = pol
	return pat

def match_ijs_lhs(pQ, t, delta_j):
	pol = Offset_list([])
	for itPIJ in pQ.Ls.ijs:
		for itSIJ in t.Ls.ijs:
			match_ijs(itPIJ, itSIJ, delta_j, pol)
	return pol

def match_ijs(itPIJ, itSIJ, delta_j, pol):
	# print itPIJ, itSIJ
	if itPIJ.i == itSIJ.i:
		if itPIJ.j + delta_j == itSIJ.j:
			pol.append(IJ(itPIJ.i, itPIJ.j))

def non_bracket_diff_length(first, second):
	length1 = 0
	length2 = 0
	i = 0
	while i < len(first):
		length1 += 1
		if first[i] == LEFTKAKKO:
			while first[i] != RIGHTKAKKO:
				i += 1
		i += 1
	i = 0
	while i < len(second):
		length2 += 1
		if second[i] == LEFTKAKKO:
			while second[i] != RIGHTKAKKO:
				i += 1
		i += 1
	return length1 - length2

##############################################################
##					Maximal modules							##
##############################################################

# lhs: Macimal_key object
# rhs: Macimal_key object
def Max_ltOp(lhs, rhs):
	if (lhs.offset_list_size < rhs.offset_list_size):
		# return True
		return 1
	elif (lhs.offset_list_size > rhs.offset_list_size):
		# return Fals
		return -1
	elif (lhs.diff_sum < rhs.diff_sum):
		# return True
		return 1
	elif (lhs.diff_sum > rhs.diff_sum):
		# return Fals
		return -1
	elif (lhs.first_sequence < rhs.first_sequence):
		# return True
		return 1
	elif (lhs.first_sequence > rhs.first_sequence):
		# return Fals
		return -1
	elif (lhs.last_sequence < rhs.last_sequence):
		# return True
		return 1
	elif (lhs.last_sequence > rhs.last_sequence):
		# return Fals
		return -1
	elif (lhs.first_diff_sum_sum < rhs.first_diff_sum_sum):
		# return True
		return 1
	elif (lhs.first_diff_sum_sum > rhs.first_diff_sum_sum):
		# return Fals
		return -1
	elif (lhs.last_diff_sum_sum < rhs.last_diff_sum_sum):
		# return True
		return 1
	elif (lhs.last_diff_sum_sum > rhs.last_diff_sum_sum):
		# return Fals
		return -1
	else:
		# return Fals
		return -1

class Maximal_key(object):
	def __init__(self):
		self.offset_list_size = 0
		self.diff_sum = 0
		self.first_sequence = 0
		self.last_sequence = 0
		self.first_diff_sum_sum = 0
		self.last_diff_sum_sum = 0

	def __eq__(self, other):
		if self.offset_list_size == other.offset_list_size:
			if self.diff_sum == other.diff_sum:
				if self.first_sequence == other.first_sequence:
					if self.last_sequence == other.last_sequence:
						if self.first_diff_sum_sum == other.first_diff_sum_sum:
							if self.last_diff_sum_sum == other.last_diff_sum_sum:
								return True
		return False

class Maximal_map_item(object):
	"""docstring for Maximal_map_item"""
	def __init__(self, Maximal_key = Maximal_key(), Pattern = Pattern()):
		self.Maximal_key = Maximal_key
		self.Pattern = Pattern

	def __eq__(self, other):
		if self.Maximal_key == other.Maximal_key and self.Pattern == other.Pattern:
			return True
		return False

class Maximal_map(object):
	def __init__(self, items = []):
		self.items = items

	def insert(self, item):
		if self.items:
			for index, elementary_pattern in enumerate(self.items):
				if Max_ltOp(elementary_pattern.Maximal_key, item.Maximal_key) > 0:
					self.items.insert(index, item)
					break
			else:
				self.items.append(item)
		else:
			self.items.append(item)

	def search(self, key):
		hit = []
		for item in self.items:
			if item.Maximal_key == key:
				hit.append(item)
		return hit

	def __iter__(self):
		return iter(self.items)

	def __getitem__(self, index):
		# print "--------"
		# print len(self.items), index
		if len(self.items) > index:
			pass
		else:
			self.items.append(Dir_s_item())
		return self.items[index]

	def __len__(self):
		return len(self.items)

def make_max_key(pattern, seqs):
	max_key = Maximal_key()
	max_key.offset_list_size = len(pattern.Ls.ijs)
	max_key.diff_sum = diff_sum(pattern.Ls.ijs, seqs)
	max_key.first_diff_sum_sum = get_global_distance(pattern.Ls.ijs[1], seqs) - get_global_distance(pattern.Ls.ijs[0], seqs)
	max_key.last_diff_sum_sum = get_global_distance(pattern.Ls.ijs[-1], seqs) - get_global_distance(pattern.Ls.ijs[-2], seqs)
	max_key.first_sequence = pattern.Ls.ijs[0].i
	max_key.last_sequence = pattern.Ls.ijs[-1].i
	return max_key

# Returns the number of brackets in a string
def count_brackets(s):
	num = 0
	for i in range(s):
		if s[i] == LEFTKAKKO:
			num+=1
	return num

def count_bracketed_length(s):
	num = 0
	i = 0
	while i < len(s):
		if s[i] == LEFTKAKKO:
			while s[i] != RIGHTKAKKO:
				num += 1
				i += 1
		i += 1
	return num

def string_search(max_motif, pat_motif):
	l1 = len(max_motif.motif) - count_bracketed_length(max_motif.motif)
	l2 = len(pat_motif.motif) - count_bracketed_length(pat_motif.motif)
	offset = max_motif.Ls.ijs[0].j - pat_motif.Ls.ijs[0].j
	for it1, it2 in zip(max_motif.Ls.ijs, pat_motif.Ls.ijs):
		if it1.i != it2.i or it1.j > it2.j or it1.j + l1 < it2.j + l2 or it1.j - it2.j != offset:
			return False
	return True

def is_subsumed(ret, pattern, maximal, max_key):
	for it in ret:
		if string_search(it.Pattern, pattern):
			return True
	return False

def add_pattern(maximal, pattern, seqs):
	max_key = make_max_key(pattern, seqs)
	ret = maximal.search(max_key)
	if not ret:
		maximal.insert(Maximal_map_item(max_key, pattern))
	elif not is_subsumed(ret, pattern, maximal, max_key):
		maximal.insert(Maximal_map_item(max_key, pattern))

def is_maximal(maximal, pattern, seqs):
	max_key = make_max_key(pattern, seqs)
	ret = maximal.search(max_key)
	# print ret
	if not ret:
		return True
	elif not is_subsumed(ret, pattern, maximal, max_key):
		return True
	return False

def get_global_distance(ij, seqs):
	sum_ = 0
	for k in range(ij.i):
		sum_ += len(seqs[k])
	return ij.j + sum_

def diff_sum(ijs, seqs):
	sum_ = 0
	len_ = len(ijs) - 1
	for k in range(len_):
		sum_ += get_global_distance(ijs[k+1], seqs) - get_global_distance(ijs[k], seqs)
	return sum_

##############################################################
##					make Dir_s modules						##
##############################################################

class Dir_s_item(object):
	def __init__(self, string = [], pattern_pairs = []):
		self.string = string
		self.pattern_pairs = pattern_pairs

	def __str__(self):
		string = "dict: "
		string += "".join(self.string) + "\n"
		for pair in self.pattern_pairs:
			string += "\t" + str(pair) + "\n"
		return string

	def delete(self, pattern):
		index, p = self.search(pattern)
		if p:
			del self.pattern_pairs[index]
		else:
			return None

	def search(self, pattern):
		for index, p in enumerate(self.pattern_pairs):
			if p == pattern:
				return index, p
		return None, None

	def is_empty(self):
		if self.pattern_pairs:
			return False
		else:
			return True

class Dir_s(object):
	"""docstring for Dir_s"""
	def __init__(self, items = []):
		self.items = items

	def insert(self, pattern):
		if self.items:
			for index, elementary_pattern in enumerate(self.items):
				if suffix_wise_less(elementary_pattern.pattern_pairs, pattern.string) > 0:
					self.items.insert(index, pattern)
					break
			else:
				self.items.append(pattern)
		else:
			self.items.append(pattern)

	def delete(self, pattern, overlap, delim):
		s = suffix(pattern.motif, overlap, delim)
		dir_index, dir_s_item = self.search(s)
		if dir_s_item and dir_index:
			dir_s_item.delete(pattern)
			if dir_s_item.is_empty():
				del self.items[dir_index]
		else:
			return None

	def search(self, suffix):
		for index, item in enumerate(self.items):
			if is_same_list(item.string, suffix):
				return index, item
		else:
			return None, Dir_s_item()

	def search_pattern(self, pattern, overlap, delim):
		s = suffix(pattern, overlap, delim)
		dir_s_index, dir_s_item = self.search(s)
		return dir_s_item.search(pattern)

	def __iter__(self):
		return iter(self.items)

	def __getitem__(self, index):
		# print "--------"
		# print len(self.items), index
		if len(self.items) > index:
			pass
		else:
			self.items.append(Dir_s_item())
		return self.items[index]

	def __len__(self):
		return len(self.items)

def suffix_wise_less(motif1, motif2):
	return 1

def suffix(pattern, overlap, delim):
	suffix = []
	count_l = 0
	i = len(pattern) - 1
	while count_l < overlap:
		if pattern[i] == LEFTKAKKO:
			while pattern[i] != RIGHTKAKKO:
				suffix.append(pattern[i])
				i -= 1
			suffix.append(pattern[i])
			count_l += 1
		elif pattern[i] != delim:
			suffix.append(pattern[i])
			count_l += 1
		else:
			suffix.append(pattern[i])
		i -= 1
	suffix = suffix[::-1]
	return suffix

def extract_s(eps, overlap_len, dir_s):
	for ep in eps:
		ep_suffix = suffix(ep.motif, overlap_len, WILDCARD)
		# print ep_suffix
		i_itDirS = 0
		# print dir_s.items
		for itDirS in dir_s:
			if is_same_list(itDirS.string, ep_suffix):
				break
			else:
				i_itDirS += 1
		if i_itDirS != len(dir_s):
			flag = False
			for i, itVector in enumerate(dir_s[i_itDirS].pattern_pairs):
				if dir_s_comp(ep.motif, itVector.motif) >= 0:
					dir_s[i_itDirS].pattern_pairs.append(ep)
					flag = True
					break
			if not flag:
				dir_s[i_itDirS].pattern_pairs.append(ep)
		else:
			pep_list = []
			pep_list.append(ep)
			dir_s.insert(Dir_s_item(ep_suffix, pep_list))
	return dir_s

def dir_s_comp(motif1, motif2):
	# print motif1, "vs", motif2
	len1 = len(motif1) - 1
	len2 = len(motif2) - 1
	length = len(motif2) if len(motif1) > len(motif2) else len(motif1)
	for i in range(length):
		if motif1[len1 -i] == WILDCARD and motif2[len2 -i] == LEFTKAKKO:
			return False
		elif motif1[len1 -i] == LEFTKAKKO and motif2[len2 -i] == WILDCARD:
			return True
		elif motif1[len1 -i] != WILDCARD and motif2[len2 -i] == WILDCARD:
			return True
		elif motif1[len1 -i] == WILDCARD and motif2[len2 -i] != WILDCARD:
			return False
		elif motif1[len1 -i] != LEFTKAKKO and motif2[len2 -i] == LEFTKAKKO:
			return False
		elif motif1[len1 -i] == LEFTKAKKO and motif2[len2 -i] != LEFTKAKKO:
			return True
	if length < len(motif1):
		return True
	elif len(motif1) < length:
		return False
	elif motif1 < motif2:
		return True
	else:
		return False

def make_dir_s(dynamic_eps, overlap_len):
	dir_s = Dir_s()
	for eps in dynamic_eps:
		dir_s = extract_s(eps, overlap_len, dir_s)
	return dir_s

##############################################################
##					make Dir_p modules						##
##############################################################

class Dir_p_item(object):
	def __init__(self, string = [], pattern_pairs = []):
		self.string = string
		self.pattern_pairs = pattern_pairs

	def __str__(self):
		return "".join(self.string)

	def delete(self, pattern):
		index, p = self.search(pattern)
		if p:
			del self.pattern_pairs[index]
		else:
			return None

	def search(self, pattern):
		for index, p in enumerate(self.pattern_pairs):
			if p == pattern:
				return index, p
		return None, None

	def is_empty(self):
		if self.pattern_pairs:
			return False
		else:
			return True

class Dir_p(object):
	"""docstring for Dir_p"""
	def __init__(self, items = []):
		self.items = items

	def insert(self, pattern):
		if self.items:
			for index, elementary_pattern in enumerate(self.items):
				if preffix_wise_less(elementary_pattern.pattern_pairs, pattern.string) > 0:
					self.items.insert(index, pattern)
					break
			else:
				self.items.append(pattern)
		else:
			self.items.append(pattern)

	def delete(self, pattern, overlap, delim):
		s = prefix(pattern.motif, overlap, delim)
		dir_index, dir_s_item = self.search(s)
		if dir_s_item and dir_index:
			dir_s_item.delete(pattern)
			if dir_s_item.is_empty():
				del self.items[dir_index]
		else:
			return None

	def search(self, prefix):
		for index, item in enumerate(self.items):
			if is_same_list(item.string, prefix):
				return index, item
		else:
			return None, Dir_p_item()

	def __iter__(self):
		return iter(self.items)

	def __getitem__(self, index):
		# print "--------"
		# print len(self.items), index
		if len(self.items) > index:
			pass
		else:
			self.items.append(Dir_p_item())
		return self.items[index]

	def __len__(self):
		return len(self.items)

def preffix_wise_less(motif1, motif2):
	return 1

def prefix(pattern, overlap, delim):
	prefix = []
	count_l = 0
	i = 0
	while count_l < overlap:
		if pattern[i] == LEFTKAKKO:
			while pattern[i] != RIGHTKAKKO:
				prefix.append(pattern[i])
				i += 1
			prefix.append(pattern[i])
			count_l += 1
		elif pattern[i] != delim:
			prefix.append(pattern[i])
			count_l += 1
		else:
			prefix.append(pattern[i])
		i += 1
	return prefix

def extract_p(eps, overlap_len, dir_p):
	for ep in eps:
		ep_prefix = prefix(ep.motif, overlap_len, WILDCARD)
		# print ep_prefix
		i_itDirP = 0
		# print dir_p.items
		for itDirP in dir_p:
			if is_same_list(itDirP.string, ep_prefix):
				break
			else:
				i_itDirP += 1
		if i_itDirP != len(dir_p):
			flag = False
			for i, itVector in enumerate(dir_p[i_itDirP].pattern_pairs):
				if dir_p_comp(ep.motif, itVector.motif) < 0:
					dir_p[i_itDirP].pattern_pairs.append(ep)
					flag = True
					break
			if not flag:
				dir_p[i_itDirP].pattern_pairs.append(ep)
		else:
			pep_list = []
			pep_list.append(ep)
			dir_p.insert(Dir_p_item(ep_prefix, pep_list))
	return dir_p

def dir_p_comp(motif1, motif2):
	# print motif1, "vs", motif2
	length = len(motif2) if len(motif1) > len(motif2) else len(motif1)
	for i in range(length):
		if motif1[i] == WILDCARD and motif2[i] == LEFTKAKKO:
			return False
		elif motif1[i] == LEFTKAKKO and motif2[i] == WILDCARD:
			return True
		elif motif1[i] != WILDCARD and motif2[i] == WILDCARD:
			return True
		elif motif1[i] == WILDCARD and motif2[i] != WILDCARD:
			return False
		elif motif1[i] != LEFTKAKKO and motif2[i] == LEFTKAKKO:
			return True
		elif motif1[i] == LEFTKAKKO and motif2[i] != LEFTKAKKO:
			return False
	flag = bracket_length_comp2(motif1, motif2)
	if flag == 1:
		return True
	elif flag == 2:
		return False
	elif (motif1 < motif2):
		return True
	else:
		return False

def make_dir_p(dynamic_eps, overlap_len):
	dir_p = Dir_p()
	for eps in dynamic_eps:
		dir_p = extract_p(eps, overlap_len, dir_p)
	return dir_p

def bracket_length_comp2(motif1, motif2):
	if (LEFTKAKKO not in motif1) and (LEFTKAKKO not in motif2):
		if len(motif1) > len(motif2):
			return 1
		elif len(motif1) < len(motif2):
			return 2
		else:
			return 0
	# bracket の処理
	return 3

def clean_up_motif(p, seqs):
	i = 0
	k = 0
	s = []
	appeared = []
	while i < len(p.motif):
		if p.motif[i] == LEFTKAKKO:
			s.append(LEFTKAKKO)
			appeared = clean_up_backet(p.Ls.ijs, k, seqs)
			for m in range(len(appeared)):
				s.append(appeared[m])
			while p.motif[i] != RIGHTKAKKO:
				i += 1
		s.append(p.motif[i])
		i += 1
		k += 1
	return s

def clean_up_soln(maximal, seqs):
	solution = []
	for ite in maximal:
		s = clean_up_motif(ite.Pattern, seqs)
		if s:
			p = Pattern(s, ite.Pattern.Ls)
			solution.append(p)
	return solution

def main():
	a = ["b", "a", ".", "d"]
	print prefix(a, 3, WILDCARD)
	print suffix(a, 2, WILDCARD)

if __name__ == "__main__":
	main()
# 
