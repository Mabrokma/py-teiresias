#!/usr/bin/env python
# coding:utf-8

from optparse import OptionParser

class Config(object):
	"""docstring for Config"""
	# optionsはoptparserのインスタンス
	def __init__(self, options):
		self.l = int(options.l)
		self.w = int(options.w)
		self.k = int(options.k)
		self.input_file = options.input_file
		self.output_file = options.output_file
		# breacketsの実装はまた後で
		self.max_brackets = -1
		self.offset_report = False
		self.convolution_length = 1
		self.max_support = -1
		self.equivalence_file = ""

def initialize_cmdline(command_line):
	parser = OptionParser()
	parser.add_option("-i", "--input", dest="input_file", help="input file", metavar="FILE")
	parser.add_option("-o", "--output", dest="output_file", default="output.txt")
	parser.add_option("-l", dest="l", default=1)
	parser.add_option("-w", dest="w", default=1)
	parser.add_option("-k", dest="k", default=1)
	command_line = command_line.split()
	(options, command_line) = parser.parse_args(command_line)
	print options
	return Config(options)

def main():
	config = initialize_cmdline("-itest.txt -oresult.txt")

if __name__ == "__main__":
	main()
