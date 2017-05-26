#!/usr/bin/env python
"""Embed signalalignment tsv into a fast5 file which is readable by nanonet

Embed most likely kmer for each event
"""
########################################################################
# File:embed_signalalign.py
#  executable: embed_signalalign.py
# Purpose: embed output from signalalign into fast5 nanopore sequencing file

# Author: Andrew Bailey
# History: 4/19/2017 Created
# NOTE:
# I will be using some ideas and code from Nanomod (https://github.com/scottgigante/nanomod)
########################################################################
from __future__ import print_function

import os
import sys
# import pysam
import csv
from timeit import default_timer as timer
from shutil import copyfile
from numpy.lib.recfunctions import append_fields, drop_fields
import numpy as np
import h5py
# from nanonet.features import *

def embed_signalalign(signalalign_tsv, events, strand, CHECK_EVERY=False, eventalign=False):
    """Embed signal align tsv file kmers into fast5 files sequencing files"""

    with open(signalalign_tsv) as tsv:
        reader = csv.reader(tsv, delimiter="\t")
        skips = 0
        stays = 0
        total = 0
        numlines = float("inf")
        if eventalign:
            next(reader)
        for line in reader:
            chromosome = np.string_(line[0])
            seq_pos = int(line[1])
            kmer = np.string_(line[15])
            name = line[3]
            template = line[4]
            event_index = int(line[5])
            prob = float(line[12])
            if total < numlines:
                if total == 0:
                    initial_event_index = event_index
                    initial_ref_index = seq_pos
                    prev_seq_pos = seq_pos
                    last_name = name
                    # total += 1
                if template == "c":
                    # break when going to complement strand
                    print("c")
                    break
                    # print(int(event_index))
                else:
                # print(line)
                # embed most probable kmer
                    if events["prob"][event_index] < prob:
                        events["prob"][event_index] = prob
                        events["kmer"][event_index] = kmer
                        if strand == "-":
                            events["seq_pos"][event_index] = initial_ref_index-seq_pos
                        else:
                            events["seq_pos"][event_index] = seq_pos-initial_ref_index
                        # TODO: fix this if needed for signalAlign
                        if seq_pos == prev_seq_pos:
                            stays += 1
                        if strand == "-":
                            if seq_pos < prev_seq_pos+1:
                                skips += 1
                        else:
                            if seq_pos > prev_seq_pos+1:
                                skips += 1
                        # this is for naomod labeling process comparison
                        if name != last_name:
                            break
                        # print (event_index, seq_pos, kmer)
                        prev_seq_pos = seq_pos
                        last_ref_index = seq_pos
                        last_name = line[3]
                        total += 1
                    if CHECK_EVERY:
                        if total > 1:
                            if int(event_index) + 1 != int(line[5]):
                                print("FAIL")
                                print(int(event_index), int(line[5]))
                                break

        # print(initial_event_index, initial_ref_index, chromosome, last_ref_index, skips, stays, total, forward)
    # print("######### EVENTS TO FEATURES ###############")
    # print(events_to_features(events)[0])
    return initial_event_index, initial_ref_index, chromosome, last_ref_index, skips, stays, total, events



def makeFast5(newPath, signalalign_tsv, fastaFile, strand, kmer=5, eventalign=False):
    """Make put in a table into a fast5 file"""
    with h5py.File(newPath, 'r+') as fast5:

        # create events
        events = fast5.get("Analyses/Basecall_1D_000/BaseCalled_template/Events").value
        # print(type(events))
        # print(events[0])
        # print(type(events))
        events = append_fields(events, 'kmer', data=['X'*kmer] * \
        events.shape[0], dtypes='<S'+str(kmer))
        events = append_fields(events, 'seq_pos', [-1] * events.shape[0])
        events = append_fields(events, 'prob', [0.0] * events.shape[0])
        events = drop_fields(events, "move")
        # print(events.shape)
        # print(events["kmer"].shape)
        # print(events["seq_pos"][1])
        initial_event_index, initial_ref_index, chromosome, last_ref_index, skips, stays, total, events = embed_signalalign(signalalign_tsv, events, strand, eventalign=eventalign)
        skipProb = skips/total
        stayProb = stays/total
        stepProb = 1-skipProb-stayProb

        b = events["seq_pos"] != -1
        events = events[b]
        events = append_fields(events, 'good_emission', events["kmer"] != 'X'*kmer)

        analysesGroup = fast5.get("Analyses")
        try:
            alignToRefGroup = analysesGroup.create_group("AlignToRef")
            eventsGroup = alignToRefGroup.create_group(("CurrentSpaceMapped_template"))
            eventsGroup.create_dataset("Events", data=events)
            # create attrs
            # dt = h5py.special_dtype(vlen=unicode)
            summaryGroup = alignToRefGroup.create_group("Summary")
            attrsGroup = summaryGroup.create_group("current_space_map_template")
            attrs = attrsGroup.attrs
            attrs.create("genome_start", initial_ref_index)
            attrs.create("genome_end", last_ref_index)
            attrs.create("sequence_length", last_ref_index - initial_ref_index)
            attrs.create("ref_name", chromosome)
            attrs.create("direction", np.string_(strand))
            attrs.create("num_skips", skips)
            attrs.create("num_stays", stays)
            attrs.create("num_events", total)
            attrs.create("skip_prob", skipProb)
            attrs.create("stay_prob", stayProb)
            attrs.create("step_prob", stepProb)

            # create alignment group
            alignGroup = analysesGroup.create_group("Alignment")
            fastaGroup = alignGroup.create_group("Aligned_template")
            fasta = open(fastaFile, "r").read()
            # print(fasta.read())
            # fasta = "asdfasdf\nATGCATGC"
            fastaGroup.create_dataset("Fasta",
                    data=np.array(fasta, dtype='<S{}'.format(len(fasta))))
            # create attrs
            summaryGroup = alignGroup.create_group("Summary")
            attrsGroup = summaryGroup.create_group("genome_mapping_template")
            attrs = attrsGroup.attrs
            attrs.create("genome", chromosome)

            print("signalAlign - Finished embedding labels into fast5 file", file=sys.stderr)

        except ValueError as e:
            print("FAST5 file already labelled. Unable to Re-Label")


def main():

    start = timer()

    signalalign_tsv = "/Users/andrewbailey/personal_git/labWork/davidHaussler/embed_signalalign/signal_align_output/tempFiles_alignment/603580b7-115c-4b57-a2f4-a3c7cddc2796_Basecall_2D_template.sm.forward.tsv"

    fast5_file = "/Users/andrewbailey/personal_git/labWork/davidHaussler/embed_signalalign/signal_align_output/tempFiles_alignment/AlexisLucattini_20160918_FNFAD24297_MN19582_sequencing_run_E_COLI_NON_MTHYLTD_R9_77950_ch215_read2042_strand1.fast5"

    eventalign_tsv = "/Users/andrewbailey/personal_git/labWork/davidHaussler/embed_signalalign/.canonical.eventalign"

    testoutput = "/Users/andrewbailey/personal_git/labWork/davidHaussler/embed_signalalign/testing/AlexisLucattini_20160918_FNFAD24297_MN19582_sequencing_run_E_COLI_NON_MTHYLTD_R9_77950_ch94_read2151_strand.fast5"

    fastaFile = "/Users/andrewbailey/personal_git/labWork/davidHaussler/embed_signalalign/signal_align_output/tempFiles_alignment/tempFiles_AlexisLucattini_20160918_FNFAD24297_MN19582_sequencing_run_E_COLI_NON_MTHYLTD_R9_77950_ch215_read2042_strand1/temp_seq_f77ad4e7-ee61-45dc-8ca5-f467bbc1ef66_Basecall_2D_template.fa"

    assert os.path.exists(signalalign_tsv)
    assert os.path.exists(fast5_file)
    assert os.path.exists(eventalign_tsv)
    # assert os.path.exists(testoutput)
    assert os.path.exists(fastaFile)

    # copyfile(fast5_file, testoutput)
    # makeFast5(fast5_file, signalalign_tsv, fastaFile, "-", eventalign=False)

    stop = timer()
    print("Running Time = {} seconds".format(stop-start), file=sys.stderr)

if __name__=="__main__":
    main()
    raise SystemExit
