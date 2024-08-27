from Bio.Blast import NCBIXML
result_handle = open('probe_GC.xml')#, 'r'
blast_records = NCBIXML.parse(result_handle)
#blast_record = blast_records.next()


#results = []
'''
with open('probe_GC.xml', 'r') as result_handle:
    blast_records = NCBIXML.parse(result_handle)
    
    for blast_record in blast_records:
        for discription in blast_record.discriptions:
            results.append({
                'title': blast_record.discriptions,
                'score': blast_record.score,
                'e-value': blast_record.e,
                'num_alignments': blast_record.num_alignments
            })
        for alignment in blast_record.alignments:
            results.append({
                'title': blast_record.discriptions,
                'score': blast_record.score,
                'e-value': blast_record.e,
                'num_alignments': blast_record.num_alignments
            })
            for hsp in alignment.hsps:
                results.append({
                    'query': blast_record.query,
                    'hit': alignment.hit_def,
                    'score': hsp.score,
                    'e_value': hsp.expect,
                    'identity': hsp.identities

                })
'''
# 打印结果
#for result in results:
#    print(result)
#print(dir(blast_record))
#print(NCBIXML.Parameters)

from Bio.Blast import NCBIXML

# 解析 BLAST 结果
with open('probe_GC.xml', 'r') as result_handle:
    blast_records = NCBIXML.parse(result_handle)

    # 获取第一个 blast_record
    blast_record = next(blast_records)  # 使用 next() 函数

    # 打印查询信息
    print(f"Query: {blast_record.query}")
    print(f"Number of alignments: {len(blast_record.alignments)}")