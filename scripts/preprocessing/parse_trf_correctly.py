import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config.paths import *

#!/usr/bin/env python3
import pandas as pd
import sys

def parse_trf_file_v2(trf_file_path: str) -> pd.DataFrame:
    """
    پارس صحیح فایل TRF با استخراج مؤتیف واقعی
    فرمت: start end period copies consensus_size percent_matches percent_indels score A B C D entropy motif sequence
    """
    
    records = []
    current_chrom = None
    
    with open(trf_file_path, 'r') as f:
        for line in f:
            line = line.strip()
            
            if line.startswith('Sequence:'):
                current_chrom = line.split(':')[1].strip()
                continue
            
            if line and line[0].isdigit():
                parts = line.split()
                
                if len(parts) >= 14:  # حداقل ۱۴ فیلد نیاز داریم
                    try:
                        # استخراج فیلدهای اصلی
                        start = int(parts[0])
                        end = int(parts[1])
                        period = float(parts[2])
                        copies = float(parts[3])
                        consensus_size = int(parts[4])
                        percent_matches = float(parts[5])
                        percent_indels = float(parts[6])
                        score = int(parts[7])
                        a = int(parts[8])
                        b = int(parts[9])
                        c = int(parts[10])
                        d = int(parts[11])
                        entropy = float(parts[12])
                        real_motif = parts[13]
                        sequence = ' '.join(parts[14:]) if len(parts) > 14 else ''
                        
                        # اگر sequence با motif شروع نشود، شاید motif بخشی از sequence است
                        if not sequence.startswith(real_motif) and len(parts) > 14:
                            # ممکن است motif تکراری در sequence باشد
                            sequence = real_motif + ' ' + ' '.join(parts[14:])
                        
                        record = {
                            'chromosome': f'chr{current_chrom}' if current_chrom and current_chrom != 'X' else 'chrX',
                            'start': start,
                            'end': end,
                            'length': end - start + 1,
                            'period': period,
                            'copies': copies,
                            'consensus_size': consensus_size,
                            'percent_matches': percent_matches,
                            'percent_indels': percent_indels,
                            'score': score,
                            'a': a,
                            'b': b,
                            'c': c,
                            'd': d,
                            'entropy': entropy,
                            'motif': real_motif,
                            'sequence': sequence
                        }
                        
                        records.append(record)
                        
                    except (ValueError, IndexError) as e:
                        print(f"⚠️ خطا در پردازش خط: {line[:50]}...")
                        continue
    
    print(f"✅ تعداد رکوردهای پارس شده: {len(records)}")
    return pd.DataFrame(records)

def validate_dataset(df):
    """اعتبارسنجی dataset"""
    print("\n🔍 اعتبارسنجی dataset:")
    print(f"  • تعداد رکوردها: {len(df):,}")
    print(f"  • طول‌های منحصربفرد: {df['length'].nunique():,}")
    print(f"  • periodهای منحصربفرد: {df['period'].nunique():,}")
    print(f"  • میانگین طول مؤتیف: {df['motif'].str.len().mean():.1f}")
    print(f"  • نمونه‌ای از مؤتیف‌ها: {', '.join(df['motif'].unique()[:5])}")
    
    # بررسی تطابق period با طول مؤتیف
    mismatched = df[abs(df['period'] - df['motif'].str.len()) > 0.5]
    if len(mismatched) > 0:
        print(f"  ⚠️  {len(mismatched)} رکورد period با طول مؤتیف تطابق ندارند")
    else:
        print("  ✅ همه periodها با طول مؤتیف تطابق دارند")

def main():
    print("🔍 شروع پارس فایل TRF (نسخه ۲)...")
    
    # پارس فایل
    df = parse_trf_file_v2('chrX.fa.2.7.7.80.10.50.500.dat')
    
    if len(df) > 0:
        # اعتبارسنجی
        validate_dataset(df)
        
        # ذخیره dataset جدید
        output_file = 'data/processed/trf_correctly_parsed_v2.tsv'
        df.to_csv(output_file, sep='\t', index=False)
        print(f"\n💾 dataset در {output_file} ذخیره شد")
        
        # نمایش نمونه
        print("\n🧪 نمونه‌ای از داده‌ها (۵ رکورد اول):")
        sample_cols = ['chromosome', 'start', 'end', 'length', 'period', 'motif', 'sequence']
        print(df[sample_cols].head().to_string())
        
        # همچنین یک فایل ساده‌تر برای تحلیل‌های بعدی
        simple_cols = ['chromosome', 'start', 'end', 'length', 'period', 'copies', 
                      'percent_matches', 'score', 'motif', 'sequence']
        simple_file = 'data/processed/trf_simple_dataset.tsv'
        df[simple_cols].to_csv(simple_file, sep='\t', index=False)
        print(f"\n📁 فایل ساده‌شده در {simple_file} ذخیره شد")
        
    else:
        print("❌ هیچ رکوردی پارس نشد!")

if __name__ == "__main__":
    main()
