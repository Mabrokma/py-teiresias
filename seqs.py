#!/usr/bin/env python
# coding:utf-8

class Sequence(object):
	"""docstring for Sequence"""
	def __init__(self, header, string):
		self.header = header
		self.string = string

	def __str__(self):
		return "".join(self.string)

	def __len__(self):
		return len(self.string)

def parse_seqs(file_name):
	f = open(file_name, "r")
	header = ""
	string = ""
	Sequence_list = []
	for line in f.readlines():
		if line.startswith(">"):
			header = line[1:]
		else:
			string = line
			Sequence_list.append(Sequence(header, [s for s in string.rstrip()]))
	f.close()
	return Sequence_list

def smallest_seq(seq_list):
	for seq in seq_list:
		print len(seq)
	return min([len(seq) for seq in seq_list])

def main():
	Sequence_list = parse_seqs("test.txt")
	for sequence in Sequence_list:
		print sequence

if __name__ == "__main__":
	main()
