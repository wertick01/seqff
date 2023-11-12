import sys
import argparse # while this is not a full-fledged product, I suggest using the program as a console utility in combination with argparse

# adding Folder_2/subfolder to the system path
sys.path.insert(0, './SeqFFLib')
# importing
from seqff import SeqFF

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--input_file', '-i', help='path to your input [.bam, .sam, .newtemp] file (for example /home/YOU/DIR/FILE.bam)', type=str, required=True)
    parser.add_argument('--output_dir', '-o', help='path to the output directory (for example /home/YOU/OUTPUT_DIR)', type=str, required=True)

    bininfo_loc="./SeqFFLib/data/pd4615-sup-0010-table2.csv"
    rdata="./SeqFFLib/data/pd4615-sup-0008-file1.rdata"

    args = parser.parse_args()

    seqff = SeqFF(
        bininfo_loc=bininfo_loc,
        input_file=args.input_file,
        rdata=rdata,
        output_loc=args.output_dir,
    )
    seqff.seqff()
    print("completed")
