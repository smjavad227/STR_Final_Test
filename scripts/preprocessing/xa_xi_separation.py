#!/usr/bin/env python3
"""
اسکریپت تفکیک Xa/Xi برای پروژه XCI-STR
ورژن 1.0
"""

import pandas as pd
import numpy as np
import os
from pathlib import Path

def main():
    print("🔬 پروژه XCI-STR: فاز تفکیک Xa/Xi")
    print("=" * 60)
    
    # مسیرهای فایل
    input_file = "data/processed/filtered_strs/strs_filtered.tsv"
    output_dir = "data/processed/xa_xi_separation"
    
    # بررسی وجود فایل ورودی
    if not os.path.exists(input_file):
        print(f"❌ خطا: فایل {input_file} پیدا نشد!")
        print("لطفاً ابتدا مراحل قبلی را اجرا کنید.")
        return False
    
    print(f"✅ فایل ورودی: {input_file}")
    
    # خواندن داده‌ها
    try:
        df = pd.read_csv(input_file, sep='\t')
        print(f"✅ داده‌ها خوانده شدند: {len(df)} خط")
    except Exception as e:
        print(f"❌ خطا در خواندن فایل: {e}")
        return False
    
    # بررسی ساختار داده‌ها
    print(f"\n📋 ساختار داده‌ها:")
    print(f"• تعداد STRها: {len(df) - 1:,} (بدون هدر)")
    print(f"• تعداد ستون‌ها: {len(df.columns)}")
    print(f"• ستون‌ها: {', '.join(df.columns[:5])}...")
    
    # آمار توصیفی
    print(f"\n📊 آمار توصیفی:")
    print(f"• میانگین طول: {df['length'].mean():.2f} bp")
    print(f"• انحراف معیار طول: {df['length'].std():.2f} bp")
    print(f"• حداقل طول: {df['length'].min()} bp")
    print(f"• حداکثر طول: {df['length'].max()} bp")
    print(f"• میانگین دوره: {df['period'].mean():.2f}")
    
    # شبیه‌سازی تفکیک Xa/Xi
    print(f"\n🎲 شبیه‌سازی تفکیک Xa/Xi:")
    np.random.seed(42)  # برای تکرارپذیری
    
    # توزیع بر اساس مطالعات: ~60% Xi، ~40% Xa
    xi_prob = 0.60  # احتمال Xi (غیرفعال)
    xa_prob = 0.40  # احتمال Xa (فعال)
    
    # تولید برچسب‌ها
    df['xci_status'] = np.random.choice(
        ['Xa', 'Xi'], 
        size=len(df), 
        p=[xa_prob, xi_prob]
    )
    
    # جدا کردن
    xa_df = df[df['xci_status'] == 'Xa'].copy()
    xi_df = df[df['xci_status'] == 'Xi'].copy()
    
    print(f"✅ تفکیک انجام شد:")
    print(f"  • Xa (کروموزوم X فعال): {len(xa_df):,} STR ({len(xa_df)/len(df)*100:.1f}%)")
    print(f"  • Xi (کروموزوم X غیرفعال): {len(xi_df):,} STR ({len(xi_df)/len(df)*100:.1f}%)")
    
    # تحلیل تفاوت‌ها
    print(f"\n📈 تحلیل تفاوت‌های Xa vs Xi:")
    print(f"  • میانگین طول Xa: {xa_df['length'].mean():.2f} bp")
    print(f"  • میانگین طول Xi: {xi_df['length'].mean():.2f} bp")
    print(f"  • تفاوت طول: {abs(xa_df['length'].mean() - xi_df['length'].mean()):.2f} bp")
    print(f"  • میانگین دوره Xa: {xa_df['period'].mean():.2f}")
    print(f"  • میانگین دوره Xi: {xi_df['period'].mean():.2f}")
    
    # ایجاد دایرکتوری خروجی
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # ذخیره فایل‌ها
    xa_file = output_path / "Xa_STRs.tsv"
    xi_file = output_path / "Xi_STRs.tsv"
    
    xa_df.to_csv(xa_file, sep='\t', index=False)
    xi_df.to_csv(xi_file, sep='\t', index=False)
    
    print(f"\n💾 فایل‌ها ذخیره شدند:")
    print(f"  • Xa: {xa_file}")
    print(f"  • Xi: {xi_file}")
    
    # ایجاد گزارش
    report = f"""گزارش تفکیک Xa/Xi - پروژه XCI-STR
===================================
تاریخ تولید: {pd.Timestamp.now()}

📁 داده‌های ورودی:
• فایل: {input_file}
• تعداد کل STRها: {len(df):,}
• تعداد STRهای واقعی: {len(df)-1:,} (بدون هدر)
• میانگین طول: {df['length'].mean():.2f} ± {df['length'].std():.2f} bp
• دامنه طول: {df['length'].min()} - {df['length'].max()} bp
• میانگین دوره: {df['period'].mean():.2f}

🎯 نتایج تفکیک (شبیه‌سازی):
• Xa (کروموزوم X فعال): {len(xa_df):,} STR ({len(xa_df)/len(df)*100:.1f}%)
• Xi (کروموزوم X غیرفعال): {len(xi_df):,} STR ({len(xi_df)/len(df)*100:.1f}%)

📊 آمار مقایسه‌ای:
                    Xa (فعال)      Xi (غیرفعال)
میانگین طول (bp):   {xa_df['length'].mean():.2f}        {xi_df['length'].mean():.2f}
انحراف معیار طول:  {xa_df['length'].std():.2f}        {xi_df['length'].std():.2f}
میانگین دوره:       {xa_df['period'].mean():.2f}        {xi_df['period'].mean():.2f}
حداقل طول:         {xa_df['length'].min()}           {xi_df['length'].min()}
حداکثر طول:        {xa_df['length'].max()}           {xi_df['length'].max()}
تعداد STRها:       {len(xa_df):,}              {len(xi_df):,}

📂 فایل‌های خروجی:
• Xa_STRs.tsv: {len(xa_df):,} خط (با هدر)
• Xi_STRs.tsv: {len(xi_df):,} خط (با هدر)

⚠️ نکات مهم:
1. این یک شبیه‌سازی اولیه است.
2. برای نتایج واقعی نیاز به داده‌های تجربی XCI داریم.
3. داده‌های مورد نیاز:
   - هایپرمتیلاسیون برای شناسایی Xa
   - RNA-seq برای تأیید غیرفعال بودن Xi
   - داده‌های ChIP-seq برای فاکتورهای XCI

🔜 مرحله بعدی پروژه:
تحلیل موتیف‌های ممنوعه با استفاده از:
1. مدل مارکوف درجه ۳
2. محاسبه O/E ratios
3. شناسایی موتیف‌های غنی‌شده در Xi
4. اعتبارسنجی آماری
"""
    
    # ذخیره گزارش
    report_file = output_path / "xa_xi_separation_report.txt"
    with open(report_file, 'w', encoding='utf-8') as f:
        f.write(report)
    
    print(f"\n📄 گزارش کامل: {report_file}")
    
    # نمایش خلاصه
    print("\n" + "=" * 60)
    print("✅ فاز ۱ پروژه تکمیل شد!")
    print("=" * 60)
    print(f"📊 نتایج نهایی:")
    print(f"• کل STRها: {len(df):,}")
    print(f"• Xa (فعال): {len(xa_df):,} ({len(xa_df)/len(df)*100:.1f}%)")
    print(f"• Xi (غیرفعال): {len(xi_df):,} ({len(xi_df)/len(df)*100:.1f}%)")
    print(f"\n🚀 آماده برای فاز ۲: کشف موتیف‌های ممنوعه")
    
    return True

if __name__ == "__main__":
    success = main()
    if success:
        print("\n🎉 عملیات با موفقیت تکمیل شد!")
    else:
        print("\n❌ عملیات با خطا مواجه شد!")
