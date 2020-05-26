#!/usr/bin/python3

import os
import subprocess
import re

name_re = re.compile(r'e(\d+)_(\d+)_(\d+)_(\d+)\.gr')

num_verts = 8
seed_set = set(range(1000))

def print_status(phase, done, total):
	print('{} {} out of {}'.format(phase, done, total), end='')

def clear_status():
	print('\r', end='');

def scan_instances():
	for ent in os.scandir('randinst/'):
		if not ent.is_file():
			continue
		match = name_re.match(ent.name)
		if not match:
			continue
		inst_verts, inst_edges = int(match.group(1)), int(match.group(2))
		if inst_verts != num_verts:
			continue
		td, seed = int(match.group(3)), int(match.group(4))
		if seed not in seed_set:
			continue
		yield (ent.path, seed, td)

done = 0
existing = set(seed for _, seed, _ in scan_instances())
for i in seed_set:
	if i in existing:
		clear_status()
		done += 1
		print_status('generating', done, len(seed_set))
		continue
	args = ['./instance_generator/main.py', '--v', str(num_verts),
				'--dir', 'randinst/', '--s', str(i)]
	rs = subprocess.call(args,
			stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
	clear_status()
	if rs != 0:
		print("seed {:5}| generator crash\n".format(i)
			+ "          | invocation: {}".format(' '.join(args)))
	done += 1
	print_status('generating', done, len(seed_set))
clear_status()

done = 0
instances = sorted(scan_instances(), key=lambda e: e[1])
for path, seed, td in instances:
	args = ['./build/td', '--solver', 'simple-pid', '--witness', path]
	proc = subprocess.run(args,
			stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
	if proc.returncode != 0:
		clear_status()
		print("seed {:5}| solver crash\n".format(seed)
			+ "          | invocation: {}".format(' '.join(args)))
		done += 1
		print_status('solving', done, len(instances))
		continue

	lines = proc.stdout.decode('ascii').rstrip().split('\n')
	if len(lines) != num_verts + 1:
		clear_status()
		print("seed {:5}| bad output format\n".format(seed)
			+ "          | invocation: {}".format(' '.join(args)))
		done += 1
		print_status('solving', done, len(instances))
		continue
	s = int(lines[0])

	if s != td:
		clear_status()
		print("seed {:5}| expected depth {}, found {}\n".format(seed, td, s)
			+ "          | invocation: {}".format(' '.join(args)))
		done += 1
		print_status('solving', done, len(instances))
		continue

	clear_status()
	done +=1
	print_status('solving', done, len(instances))
clear_status()
