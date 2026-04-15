#!/usr/bin/env python3
"""
Script para gerar um HMM proteico válido para testes do pipeline not2code.
Requer hmmbuild (pacote HMMER) instalado no ambiente.

Uso:
    micromamba activate nextflow
    python3 assets/generate_test_hmm.py

Ou diretamente com singularity:
    singularity exec \
        /home/gmdazevedo/scratch/singularity_images/depot.galaxyproject.org-singularity-hmmer-3.4--hdbdd923_1.img \
        hmmbuild --amino assets/test_protein.hmm /dev/stdin << 'EOF'
# STOCKHOLM 1.0
seq1  MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAAVRESVPSLL
//
EOF
"""

import subprocess
import sys
import os
import textwrap

# Sequência proteica de teste (Globina humana - proteína bem caracterizada)
# Suficientemente longa para TransDecoder e para hmmbuild criar um HMM válido
TEST_SEQUENCE = """\
MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLD
NLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH
""".replace('\n', '')

STOCKHOLM_CONTENT = f"""# STOCKHOLM 1.0
#=GF ID   test_protein_hmm
#=GF AC   TEST0001.1
#=GF DE   Test protein HMM for not2code pipeline - amino acid profile
test_seq  {TEST_SEQUENCE}
//
"""

def main():
    # Caminho para o arquivo de saída
    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_hmm = os.path.join(script_dir, "test_protein.hmm")
    sto_file = os.path.join(script_dir, "test_protein.sto")

    print(f"Gerando HMM proteico de teste em: {output_hmm}")

    # Escrever o alinhamento Stockholm
    with open(sto_file, 'w') as f:
        f.write(STOCKHOLM_CONTENT)
    print(f"Alinhamento Stockholm criado: {sto_file}")

    # Executar hmmbuild
    cmd = ["hmmbuild", "--amino", output_hmm, sto_file]
    print(f"Executando: {' '.join(cmd)}")
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"ERRO ao executar hmmbuild:")
        print(result.stderr)
        sys.exit(1)
    
    print(f"HMM proteico gerado com sucesso: {output_hmm}")
    
    # Verificar o arquivo gerado
    if os.path.exists(output_hmm):
        size = os.path.getsize(output_hmm)
        print(f"Tamanho do arquivo: {size} bytes")
        
        # Mostrar o cabeçalho do HMM
        with open(output_hmm) as f:
            lines = f.readlines()[:8]
            print("Primeiras linhas do HMM:")
            for line in lines:
                print(f"  {line}", end='')
    
    # Limpar arquivo temporário
    os.remove(sto_file)
    print("\nConcluído! Use este arquivo em test.config como pfam_db.")
    print(f'pfam_db = "{output_hmm}"')

if __name__ == "__main__":
    main()
