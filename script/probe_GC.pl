#!/usr/bin/perl

use strict;
use warnings;
use Bio::SearchIO;

# 指定要解析的 XML 檔案
my $blast_xml = "probe_GC.xml";

# 創建 Bio::SearchIO 對象
my $searchio = Bio::SearchIO->new(
    -format => 'blastxml',
    -file   => $blast_xml
);

# 打開輸出CSV檔案
open my $out_fh, '>', 'probe_GC_parsed.csv' or die "Cannot open output file: $!";

# 寫入CSV標題行
print $out_fh "Query_ID,Hit_ID,E-value,Score,Identity,Alignment_Length\n";

# 迭代每個查詢結果
while (my $result = $searchio->next_result) {
    # 迭代每個命中
    while (my $hit = $result->next_hit) {
        # 迭代每個 HSP
        while (my $hsp = $hit->next_hsp) {
            # 將結果寫入CSV格式
            print $out_fh join(",",
                $result->query_name,
                $hit->name,
                $hsp->evalue,
                $hsp->score,
                $hsp->percent_identity,
                $hsp->length("total")
            ), "\n";
        }
    }
}

close $out_fh;