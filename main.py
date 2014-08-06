#!/usr/bin/env python
# coding:utf-8

from config import *
from seqs import *
from pattern import *
from convolution import *

import copy

import sys

LITERAL = 1

def main():
	command_line = sys.argv[1:]
	config = initialize_cmdline(command_line)
	seqs = parse_seqs(config.input_file)
	results = teiresias(config, seqs)
	for result in results:
		print "%s\t%s\t%s" % (result[0], result[1], result[2])

def teiresias(config, seqs):
	print config.l, config.w, config.k
	if config.l <= config.w and config.w <= smallest_seq(seqs) and config.l >= 2:
		eps = scan(seqs, config)
		print "--------------"
		for ep in eps:
			print ep
		results = convolute(eps, config, seqs)
		print "--------------"
		title = """##########################################################
#                                                        #
#                       FINAL RESULTS                    #
#                                                        #
##########################################################
"""
		f = open(config.output_file, "w")
		f.write(title)
		for result in results:
			f.write("%s\t%s\t%s\n" % (result[0], result[1], result[2]))
		f.close()
		return results
	else:
		print "error occured. check parameters."

def scan(seqs, config):
	print "Scanning..."

	subst_alpha = []
	# equivalence file handled by alphabet.h
	bits = perm_bits(config.w - 1, config.l - 1)
	bits = ["1" + bit for bit in bits]

	eps = Elementary_patterns()
	final_eps = Elementary_patterns()
	for current_bit_mask in bits:
		current_bit_mask = trim0s(current_bit_mask)

		for i in range(len(seqs)):
			# print "i", i
			end = len(seqs[i].string) - len(current_bit_mask) + 1
			for j in range(end):
				# print "j", j
				eps = make_and_insert(subst_alpha, eps, current_bit_mask, seqs, i, j)

		# print "-----------------"

		eps.supported(config)

		for eps_it in eps:
			final_eps.insert(eps_it)
	return eps

class LoopException(Exception):
	pass

class DebugException(Exception):
	pass

def convolute(eps, config, seqs):
	dynamic_eps2 = []
	dynamic_eps2.append(eps)

	print "Convoluting..."
	print "==================================="
	OVERLAP_LEN = config.l - LITERAL
	dir_p = make_dir_p(dynamic_eps2, OVERLAP_LEN)
	dir_s = make_dir_s(dynamic_eps2, OVERLAP_LEN)

	maximal = Maximal_map()
	pattern_stack = []
	w = []
	smt = []
	u_pepvec = []
	count = 0
	for eps in dynamic_eps2:
		for ep in sorted(eps, cmp = lambda x, y:alphabetical(y.motif, x.motif)):
			p = copy.deepcopy(ep)
			print "--------------"
			print "Push1 %s" % p
			pattern_stack.append(p)
			flag = True
			# label start
			while flag:
				try:
					# print "============ try ==========="
					if len(pattern_stack) == 0:
						print "remove1", p
						remove_entries(p, dir_p, dir_s, OVERLAP_LEN, WILDCARD)
						# break
						flag = False
						raise LoopException
					t = pattern_stack[-1]

					w = prefix(t.motif, OVERLAP_LEN, WILDCARD)
					index, u_pepvec = dir_s.search(w)
					# print u_pepvec
					for itU in u_pepvec.pattern_pairs:
						pQ = itU
						r = left_convolute(pQ, t, w)
						# print r
						# print "r offset : " + str(r.Ls)
						# print "pQoffset : " + str(pQ.Ls)
						# print "t offset : " + str(t.Ls)
						if config.equivalence_file == "" and len(r.Ls.ijs) == len(t.Ls.ijs):
							if pattern_stack:
								poped = pattern_stack[-1]
								print "Pop1 %s" % poped
								pattern_stack.pop()
						# print "suffix_is_maximal"
						if len(r.Ls.ijs) >= config.k and is_maximal(maximal, r, seqs):
							# print "True"
							print "Push2 %s" % r
							pattern_stack.append(r)
							# goto start
							flag = True
							raise LoopException
						else:
							# print "False"
							pass

					w = suffix(t.motif, OVERLAP_LEN, WILDCARD)
					index, u_pepvec = dir_p.search(w)
					# print w, u_pepvec
					for itU in u_pepvec.pattern_pairs:
						pQ = itU
						r = right_convolute(t, pQ, w)
						# print r
						# print "r offset : " + str(r.Ls)
						# print "pQoffset : " + str(pQ.Ls)
						# print "t offset : " + str(t.Ls)
						if config.equivalence_file == "" and len(r.Ls.ijs) == len(t.Ls.ijs):
							if pattern_stack:
								poped = pattern_stack[-1]
								print "Pop2 %s" % poped
								pattern_stack.pop()
						# print "prefix_is_maximal?"
						if len(r.Ls.ijs) >= config.k and is_maximal(maximal, r, seqs):
							# print "True"
							print "Push3 %s" % r
							pattern_stack.append(r)
							# goto start
							flag = True
							raise LoopException
						else:
							pass
							# print "False"

					if not pattern_stack:
						print "remove2", t
						remove_entries(t, dir_p, dir_s, OVERLAP_LEN, WILDCARD)
						print "touru"
						# break
						flag = False
						raise LoopException
					add_pattern(maximal, t, seqs)
					poped = pattern_stack[-1]
					pattern_stack.pop()
					print "Pop3 %s" % poped
				except LoopException:
					pass
	# ここらへんにちょっとbracksの処理あり
	print "-----------------------"
	soln_vec = clean_up_soln(maximal, seqs)
	results = []
	for sol in soln_vec:
		len_of_offset = len(sol.Ls)
		len_of_cover = len(sol.Ls.cover_strings())
		motif = str(sol)
		results.append((len_of_offset, len_of_cover, motif))
	return results

if __name__ == "__main__":
	main()
