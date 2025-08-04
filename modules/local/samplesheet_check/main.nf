#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process SAMPLESHEET_CHECK {
    tag "$samplesheet"
    label 'process_single'
    
    container 'quay.io/biocontainers/python:3.8.3'
    
    input:
    path samplesheet
    
    output:
    path '*.csv', emit: csv
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    """
    #!/usr/bin/env python3
    
    import csv
    import logging
    import sys
    import os
    from pathlib import Path
    
    def check_samplesheet(file_in, file_out):

        sample_mapping_dict = {}
        
        # Verificar se o arquivo de entrada existe
        if not os.path.exists(file_in):
            print(f"ERROR: Arquivo samplesheet não encontrado: {file_in}")
            sys.exit(1)
        
        # Verificar se o arquivo não está vazio
        if os.path.getsize(file_in) == 0:
            print(f"ERROR: Arquivo samplesheet está vazio: {file_in}")
            sys.exit(1)
        
        print(f"INFO: Validando samplesheet: {file_in}")
        
        with open(file_in, "r") as fin:
            # Verificar se o arquivo tem pelo menos uma linha
            first_line = fin.readline().strip()
            if not first_line:
                print("ERROR: Samplesheet está vazio ou só contém linhas em branco")
                sys.exit(1)
            
            # Voltar ao início do arquivo
            fin.seek(0)
            
            try:
                reader = csv.DictReader(fin)
            except Exception as e:
                print(f"ERROR: Erro ao ler CSV: {e}")
                sys.exit(1)
            
            # Verificar cabeçalhos obrigatórios
            required_columns = ["sample", "gtf"]
            if not reader.fieldnames:
                print("ERROR: Nenhum cabeçalho encontrado no arquivo CSV")
                sys.exit(1)
                
            missing_columns = [col for col in required_columns if col not in reader.fieldnames]
            if missing_columns:
                print(f"ERROR: Colunas obrigatórias não encontradas: {missing_columns}")
                print(f"Colunas encontradas: {list(reader.fieldnames)}")
                sys.exit(1)
            
            print(f"INFO: Cabeçalhos válidos encontrados: {list(reader.fieldnames)}")
            
            # Processar linhas
            line_count = 0
            for line_num, row in enumerate(reader, start=2):
                sample = row["sample"].strip()
                gtf = row["gtf"].strip()
                
                # Verificar se campos não estão vazios
                if not sample:
                    print(f"ERROR: Linha {line_num}: Campo 'sample' está vazio")
                    sys.exit(1)
                if not gtf:
                    print(f"ERROR: Linha {line_num}: Campo 'gtf' está vazio")
                    sys.exit(1)
                
                # Verificar caracteres inválidos no nome da amostra
                invalid_chars = [' ', '/', '\\\\', ':', '*', '?', '"', '<', '>', '|']
                if any(char in sample for char in invalid_chars):
                    print(f"ERROR: Linha {line_num}: Nome da amostra '{sample}' contém caracteres inválidos")
                    print(f"Caracteres não permitidos: {invalid_chars}")
                    sys.exit(1)
                
                # Verificar duplicatas
                if sample in sample_mapping_dict:
                    print(f"ERROR: Linha {line_num}: Sample '{sample}' duplicado")
                    print(f"Primeira ocorrência na linha com GTF: {sample_mapping_dict[sample]}")
                    sys.exit(1)
                
                # Verificar se o arquivo GTF tem extensão válida
                if not gtf.lower().endswith(('.gtf', '.gff', '.gff3')):
                    print(f"WARNING: Linha {line_num}: Arquivo '{gtf}' não tem extensão GTF/GFF esperada")
                
                sample_mapping_dict[sample] = gtf
                line_count += 1
            
            if line_count == 0:
                print("ERROR: Nenhuma amostra válida encontrada na samplesheet")
                sys.exit(1)
            
            print(f"INFO: {line_count} amostras válidas encontradas")
        
        # Escrever samplesheet validada
        print(f"INFO: Escrevendo samplesheet validada: {file_out}")
        with open(file_out, "w", newline='') as fout:
            writer = csv.writer(fout)
            writer.writerow(["sample", "gtf"])
            for sample, gtf in sample_mapping_dict.items():
                writer.writerow([sample, gtf])
        
        print("INFO: Validação da samplesheet concluída com sucesso!")
        
        # Imprimir resumo
        print("\\n=== RESUMO ===")
        print(f"Total de amostras: {len(sample_mapping_dict)}")
        print("Amostras encontradas:")
        for i, (sample, gtf) in enumerate(sample_mapping_dict.items(), 1):
            print(f"  {i:2d}. {sample} -> {gtf}")
    
    # Executar validação
    check_samplesheet("${samplesheet}", "samplesheet.valid.csv")
    
    # Criar arquivo de versões
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write('    python: "3.8.3"\\n')
    """
    
    stub:
    """
    echo "sample,gtf" > samplesheet.valid.csv
    echo "test_sample,test.gtf" >> samplesheet.valid.csv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "3.8.3"
    END_VERSIONS
    """
}