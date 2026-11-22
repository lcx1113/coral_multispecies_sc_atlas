import re

def parse_gtf_attributes(attr_str):
    attrs = {}
    for k,v in re.findall(r'(\S+)\s+"([^"]+)"', attr_str):
        attrs[k] = v
    return attrs

def format_gtf_line(fields, attrs):
    attr_field = ' '.join([f'{k} "{v}";' for k,v in attrs.items()])
    return '\t'.join(fields[:8]) + '\t' + attr_field

def merge_intervals(intervals):
    """ 合并重叠区间 """
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = [intervals[0]]
    for start, end in intervals[1:]:
        last_start, last_end = merged[-1]
        if start <= last_end+1:
            merged[-1][1] = max(last_end, end)
        else:
            merged.append([start,end])
    return merged

def load_gtf_blocks(gtf_file):
    # transcript: {(gene_id,transcript_id): [fields,attrs]} # 只保留transcript条
    # exon: {(gene_id,transcript_id): [[start,end,fields,attrs], ...]}
    transcripts = {}
    exons = {}
    with open(gtf_file) as f:
        for line in f:
            if line.startswith('#') or not line.strip(): continue
            fields = line.strip().split('\t')
            feature = fields[2]
            attrs = parse_gtf_attributes(fields[8])
            gene_id = attrs.get('gene_id')
            transcript_id = attrs.get('transcript_id')
            # gene变成transcript！
            if feature in ('gene','transcript'):
                fields2 = fields[:]
                fields2[2] = 'transcript'
                key = (gene_id, transcript_id)
                transcripts[key] = (fields2, attrs)
            if feature=='exon':
                key = (gene_id, transcript_id)
                if key not in exons:
                    exons[key] = []
                exons[key].append([int(fields[3]),int(fields[4]), fields, attrs])
    return transcripts, exons

def main(extend_gtf, orig_gtf, output_gtf):
    ## extend_gtf: 获取主区段-transcript
    transcripts, _ = load_gtf_blocks(extend_gtf)
    ## orig_gtf: 获取所有exon
    _, orig_exons = load_gtf_blocks(orig_gtf)
    gtf_lines = []

    for (gene_id, tx_id) in transcripts:
        # 1. transcript输出
        fields, attrs = transcripts[(gene_id, tx_id)]
        gtf_lines.append(format_gtf_line(fields, attrs))

        # 2. 搜索orig中所有exon
        exon_blocks = []
        if (gene_id, tx_id) in orig_exons:
            for s,e,fields_e,attrs_e in sorted(orig_exons[(gene_id, tx_id)]):
                exon_blocks.append([s,e,fields_e,attrs_e])
                gtf_lines.append(format_gtf_line(fields_e,attrs_e))

        # orphan: orig无exon,则输出全区间为一个exon
        if not exon_blocks:
            t_start = int(fields[3])
            t_end = int(fields[4])
            exon_fields = fields[:]
            exon_fields[2] = "exon"
            exon_attrs = attrs.copy()
            exon_attrs["exon_id"] = "orphan_exon"
            exon_fields[3] = str(t_start)
            exon_fields[4] = str(t_end)
            gtf_lines.append(format_gtf_line(exon_fields,exon_attrs))
            continue

        # 检查需不需要补3'端UTR为exon
        # 找末端exon的位置，仅补最末端到transcript末端
        t_start = int(fields[3])
        t_end = int(fields[4])
        strand = fields[6]

        # 找最末exon
        if strand == "+":
            exon_blocks_sorted = sorted(exon_blocks, key=lambda x: x[1], reverse=True)
            last_exon = exon_blocks_sorted[0]
            last_exon_end = last_exon[1]
            if last_exon_end < t_end:
                # 必须有区间，否则不补
                new_exon_start = last_exon_end + 1
                new_exon_end = t_end
                if new_exon_start <= new_exon_end:
                    exon_fields = fields[:]
                    exon_fields[2] = "exon"
                    exon_fields[3] = str(new_exon_start)
                    exon_fields[4] = str(new_exon_end)
                    exon_attrs = attrs.copy()
                    exon_attrs["exon_id"] = "extended_3putr"
                    gtf_lines.append(format_gtf_line(exon_fields, exon_attrs))
        elif strand == "-":
            exon_blocks_sorted = sorted(exon_blocks, key=lambda x: x[0])
            last_exon = exon_blocks_sorted[0]
            last_exon_start = last_exon[0]
            if t_start < last_exon_start:
                new_exon_start = t_start
                new_exon_end = last_exon_start - 1
                if new_exon_start <= new_exon_end:
                    exon_fields = fields[:]
                    exon_fields[2] = "exon"
                    exon_fields[3] = str(new_exon_start)
                    exon_fields[4] = str(new_exon_end)
                    exon_attrs = attrs.copy()
                    exon_attrs["exon_id"] = "extended_3putr"
                    gtf_lines.append(format_gtf_line(exon_fields, exon_attrs))
    # 输出
    with open(output_gtf,"w") as f:
        for line in gtf_lines:
            f.write(line+'\n')

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--extend_gtf", required=True, help="GTF: extended gene as transcript")
    parser.add_argument("--orig_gtf", required=True, help="GTF: original exons")
    parser.add_argument("--output", required=True, help="Output merged GTF")
    args = parser.parse_args()
    main(args.extend_gtf, args.orig_gtf, args.output)