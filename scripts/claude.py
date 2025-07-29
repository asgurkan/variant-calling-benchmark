import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path
import os

# Türkçe font ve stil ayarları
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['axes.facecolor'] = 'white'

class MetricsVisualizer:
    def __init__(self, input_folder="input", output_folder="output"):
        self.input_folder = Path(input_folder)
        self.output_folder = Path(output_folder)
        self.output_folder.mkdir(exist_ok=True)
        
        # Renk paleti - modern ve estetik
        self.color_palette = {
            'Medaka': '#2E86AB',      # Mavi
            'Longshot': '#A23B72',    # Mor
            'good': '#43AA8B',        # Yeşil
            'medium': '#F18F01',      # Turuncu
            'bad': '#C73E1D'          # Kırmızı
        }
        
        # Stil ayarları
        sns.set_style("whitegrid", {
            'axes.grid': True,
            'axes.edgecolor': '.8',
            'grid.color': '.92'
        })
        
    def load_data(self, filename):
        """CSV dosyasını yükle"""
        filepath = self.input_folder / filename
        try:
            df = pd.read_csv(filepath)
            print(f"✅ {filename} başarıyla yüklendi - {len(df)} satır")
            return df
        except Exception as e:
            print(f"❌ Hata: {filename} yüklenemedi - {e}")
            return None
    
    def create_performance_comparison(self, df):
        """Araçlar arası performans karşılaştırması"""
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        fig.suptitle('🔬 Araç Performans Karşılaştırması', fontsize=20, fontweight='bold', y=0.98)
        
        # Precision karşılaştırması
        ax1 = axes[0, 0]
        precision_data = df.pivot_table(values='Precision', index='Quality', columns='Tool', aggfunc='mean')
        sns.heatmap(precision_data, annot=True, fmt='.4f', cmap='RdYlGn', 
                   ax=ax1, cbar_kws={'label': 'Precision'})
        ax1.set_title('📊 Precision Değerleri', fontsize=14, fontweight='bold')
        ax1.set_xlabel('Araç', fontweight='bold')
        ax1.set_ylabel('Kalite', fontweight='bold')
        
        # Recall karşılaştırması
        ax2 = axes[0, 1]
        recall_data = df.pivot_table(values='Recall', index='Quality', columns='Tool', aggfunc='mean')
        sns.heatmap(recall_data, annot=True, fmt='.4f', cmap='RdYlGn', 
                   ax=ax2, cbar_kws={'label': 'Recall'})
        ax2.set_title('📈 Recall Değerleri', fontsize=14, fontweight='bold')
        ax2.set_xlabel('Araç', fontweight='bold')
        ax2.set_ylabel('Kalite', fontweight='bold')
        
        # F1-Score karşılaştırması
        ax3 = axes[1, 0]
        f1_data = df.pivot_table(values='F1-score', index='Quality', columns='Tool', aggfunc='mean')
        sns.heatmap(f1_data, annot=True, fmt='.4f', cmap='RdYlGn', 
                   ax=ax3, cbar_kws={'label': 'F1-Score'})
        ax3.set_title('🎯 F1-Score Değerleri', fontsize=14, fontweight='bold')
        ax3.set_xlabel('Araç', fontweight='bold')
        ax3.set_ylabel('Kalite', fontweight='bold')
        
        # Genel performans radarı
        ax4 = axes[1, 1]
        tools = df['Tool'].unique()
        metrics = ['Precision', 'Recall', 'F1-score']
        
        for tool in tools:
            tool_data = df[df['Tool'] == tool][metrics].mean()
            ax4.plot(metrics, tool_data.values, 'o-', linewidth=3, markersize=8, 
                    label=tool, color=self.color_palette[tool])
            ax4.fill_between(metrics, tool_data.values, alpha=0.2, 
                           color=self.color_palette[tool])
        
        ax4.set_title('🎭 Genel Performans Karşılaştırması', fontsize=14, fontweight='bold')
        ax4.set_ylabel('Skor', fontweight='bold')
        ax4.legend(frameon=True, fancybox=True, shadow=True)
        ax4.grid(True, alpha=0.3)
        ax4.set_ylim(0, 1.1)
        
        plt.tight_layout()
        return fig
    
    def create_quality_analysis(self, df):
        """Kalite bazlı analiz"""
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        fig.suptitle('🌟 Kalite Bazlı Detaylı Analiz', fontsize=20, fontweight='bold', y=0.98)
        
        # Kalite dağılımı
        ax1 = axes[0, 0]
        quality_counts = df['Quality'].value_counts()
        colors = [self.color_palette[q] for q in quality_counts.index]
        wedges, texts, autotexts = ax1.pie(quality_counts.values, labels=quality_counts.index, 
                                          autopct='%1.1f%%', startangle=90, colors=colors,
                                          textprops={'fontsize': 12, 'fontweight': 'bold'})
        ax1.set_title('🥧 Kalite Dağılımı', fontsize=14, fontweight='bold')
        
        # Kaliteye göre F1-Score dağılımı
        ax2 = axes[0, 1]
        sns.boxplot(data=df, x='Quality', y='F1-score', hue='Tool', ax=ax2,
                   palette=[self.color_palette[tool] for tool in df['Tool'].unique()])
        ax2.set_title('📦 Kaliteye Göre F1-Score Dağılımı', fontsize=14, fontweight='bold')
        ax2.set_xlabel('Kalite', fontweight='bold')
        ax2.set_ylabel('F1-Score', fontweight='bold')
        ax2.legend(title='Araç', frameon=True, fancybox=True, shadow=True)
        
        # Hata analizi (FP ve FN)
        ax3 = axes[1, 0]
        error_data = df[['Tool', 'Quality', 'FP', 'FN']].melt(
            id_vars=['Tool', 'Quality'], var_name='Error_Type', value_name='Count')
        sns.barplot(data=error_data, x='Quality', y='Count', hue='Error_Type', ax=ax3,
                   palette={'FP': '#FF6B6B', 'FN': '#4ECDC4'})
        ax3.set_title('⚠️ Hata Analizi (False Positive & False Negative)', fontsize=14, fontweight='bold')
        ax3.set_xlabel('Kalite', fontweight='bold')
        ax3.set_ylabel('Hata Sayısı', fontweight='bold')
        ax3.legend(title='Hata Tipi', frameon=True, fancybox=True, shadow=True)
        ax3.set_yscale('log')  # Log scale for better visibility
        
        # Precision vs Recall scatter
        ax4 = axes[1, 1]
        for tool in df['Tool'].unique():
            tool_data = df[df['Tool'] == tool]
            for quality in tool_data['Quality'].unique():
                quality_data = tool_data[tool_data['Quality'] == quality]
                ax4.scatter(quality_data['Precision'], quality_data['Recall'], 
                          s=200, alpha=0.7, color=self.color_palette[tool],
                          marker='o' if tool == 'Medaka' else 's',
                          label=f'{tool}-{quality}' if quality_data.index[0] == df.index[0] else "")
        
        ax4.set_title('🎯 Precision vs Recall Analizi', fontsize=14, fontweight='bold')
        ax4.set_xlabel('Precision', fontweight='bold')
        ax4.set_ylabel('Recall', fontweight='bold')
        ax4.legend(bbox_to_anchor=(1.05, 1), loc='upper left', frameon=True, fancybox=True, shadow=True)
        ax4.grid(True, alpha=0.3)
        ax4.set_xlim(0, 1.1)
        ax4.set_ylim(0, 1.1)
        
        plt.tight_layout()
        return fig
    
    def create_detailed_metrics_table(self, df):
        """Detaylı metrik tablosu"""
        fig, ax = plt.subplots(figsize=(14, 8))
        fig.suptitle('📋 Detaylı Metrik Tablosu', fontsize=18, fontweight='bold', y=0.95)
        
        # Tabloyu hazırla
        table_data = df.round(4)
        
        # Renk kodlaması için normalizasyon
        def get_color_for_value(value, metric):
            if metric in ['Precision', 'Recall', 'F1-score']:
                if value >= 0.8:
                    return '#43AA8B'  # Yeşil
                elif value >= 0.5:
                    return '#F18F01'  # Turuncu
                else:
                    return '#C73E1D'  # Kırmızı
            return 'white'
        
        # Tabloyu oluştur
        table = ax.table(cellText=table_data.values,
                        colLabels=table_data.columns,
                        cellLoc='center',
                        loc='center',
                        bbox=[0, 0, 1, 1])
        
        table.auto_set_font_size(False)
        table.set_fontsize(10)
        table.scale(1, 2)
        
        # Hücre renklerini ayarla
        for i in range(len(table_data)):
            for j, col in enumerate(table_data.columns):
                if col in ['Precision', 'Recall', 'F1-score']:
                    color = get_color_for_value(table_data.iloc[i, j], col)
                    table[(i+1, j)].set_facecolor(color)
                    table[(i+1, j)].set_text_props(weight='bold', color='white' if color != 'white' else 'black')
        
        # Başlık satırını stillendir
        for j in range(len(table_data.columns)):
            table[(0, j)].set_facecolor('#2E86AB')
            table[(0, j)].set_text_props(weight='bold', color='white')
        
        ax.axis('off')
        return fig
    
    def create_summary_dashboard(self, df):
        """Özet dashboard"""
        fig = plt.figure(figsize=(20, 12))
        gs = fig.add_gridspec(3, 4, hspace=0.3, wspace=0.3)
        fig.suptitle('🚀 Kapsamlı Metrik Dashboard', fontsize=24, fontweight='bold', y=0.95)
        
        # En iyi performansları bul
        best_precision = df.loc[df['Precision'].idxmax()]
        best_recall = df.loc[df['Recall'].idxmax()]
        best_f1 = df.loc[df['F1-score'].idxmax()]
        
        # KPI kartları
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.text(0.5, 0.7, f"🏆 En İyi Precision", ha='center', va='center', 
                fontsize=14, fontweight='bold', transform=ax1.transAxes)
        ax1.text(0.5, 0.5, f"{best_precision['Precision']:.4f}", ha='center', va='center', 
                fontsize=20, fontweight='bold', color=self.color_palette[best_precision['Tool']], 
                transform=ax1.transAxes)
        ax1.text(0.5, 0.3, f"{best_precision['Tool']} - {best_precision['Quality']}", 
                ha='center', va='center', fontsize=12, transform=ax1.transAxes)
        ax1.set_xlim(0, 1)
        ax1.set_ylim(0, 1)
        ax1.axis('off')
        ax1.add_patch(plt.Rectangle((0.05, 0.05), 0.9, 0.9, fill=False, edgecolor='gray', linewidth=2))
        
        ax2 = fig.add_subplot(gs[0, 1])
        ax2.text(0.5, 0.7, f"🎯 En İyi Recall", ha='center', va='center', 
                fontsize=14, fontweight='bold', transform=ax2.transAxes)
        ax2.text(0.5, 0.5, f"{best_recall['Recall']:.4f}", ha='center', va='center', 
                fontsize=20, fontweight='bold', color=self.color_palette[best_recall['Tool']], 
                transform=ax2.transAxes)
        ax2.text(0.5, 0.3, f"{best_recall['Tool']} - {best_recall['Quality']}", 
                ha='center', va='center', fontsize=12, transform=ax2.transAxes)
        ax2.set_xlim(0, 1)
        ax2.set_ylim(0, 1)
        ax2.axis('off')
        ax2.add_patch(plt.Rectangle((0.05, 0.05), 0.9, 0.9, fill=False, edgecolor='gray', linewidth=2))
        
        ax3 = fig.add_subplot(gs[0, 2])
        ax3.text(0.5, 0.7, f"⭐ En İyi F1-Score", ha='center', va='center', 
                fontsize=14, fontweight='bold', transform=ax3.transAxes)
        ax3.text(0.5, 0.5, f"{best_f1['F1-score']:.4f}", ha='center', va='center', 
                fontsize=20, fontweight='bold', color=self.color_palette[best_f1['Tool']], 
                transform=ax3.transAxes)
        ax3.text(0.5, 0.3, f"{best_f1['Tool']} - {best_f1['Quality']}", 
                ha='center', va='center', fontsize=12, transform=ax3.transAxes)
        ax3.set_xlim(0, 1)
        ax3.set_ylim(0, 1)
        ax3.axis('off')
        ax3.add_patch(plt.Rectangle((0.05, 0.05), 0.9, 0.9, fill=False, edgecolor='gray', linewidth=2))
        
        # Genel istatistikler
        ax4 = fig.add_subplot(gs[0, 3])
        total_tp = df['TP'].sum()
        total_fp = df['FP'].sum()
        total_fn = df['FN'].sum()
        ax4.text(0.5, 0.8, f"📊 Genel İstatistikler", ha='center', va='center', 
                fontsize=14, fontweight='bold', transform=ax4.transAxes)
        ax4.text(0.5, 0.6, f"TP: {total_tp:,.0f}", ha='center', va='center', 
                fontsize=12, color='green', transform=ax4.transAxes)
        ax4.text(0.5, 0.45, f"FP: {total_fp:,.0f}", ha='center', va='center', 
                fontsize=12, color='red', transform=ax4.transAxes)
        ax4.text(0.5, 0.3, f"FN: {total_fn:,.0f}", ha='center', va='center', 
                fontsize=12, color='orange', transform=ax4.transAxes)
        ax4.set_xlim(0, 1)
        ax4.set_ylim(0, 1)
        ax4.axis('off')
        ax4.add_patch(plt.Rectangle((0.05, 0.05), 0.9, 0.9, fill=False, edgecolor='gray', linewidth=2))
        
        # Ana grafikler
        ax5 = fig.add_subplot(gs[1, :2])
        tool_performance = df.groupby(['Tool', 'Quality'])[['Precision', 'Recall', 'F1-score']].mean()
        tool_performance.plot(kind='bar', ax=ax5, color=['#FF6B6B', '#4ECDC4', '#45B7D1'], 
                             width=0.8, alpha=0.8)
        ax5.set_title('📈 Araç ve Kalite Bazlı Performans', fontsize=14, fontweight='bold')
        ax5.legend(frameon=True, fancybox=True, shadow=True)
        ax5.tick_params(axis='x', rotation=45)
        
        ax6 = fig.add_subplot(gs[1, 2:])
        correlation_matrix = df[['TP', 'FP', 'FN', 'Precision', 'Recall', 'F1-score']].corr()
        sns.heatmap(correlation_matrix, annot=True, cmap='RdBu_r', center=0, ax=ax6,
                   square=True, fmt='.2f')
        ax6.set_title('🔗 Metrik Korelasyon Matrisi', fontsize=14, fontweight='bold')
        
        # Alt panel - trend analizi
        ax7 = fig.add_subplot(gs[2, :])
        df_sorted = df.sort_values(['Quality', 'Tool'])
        x_pos = range(len(df_sorted))
        
        ax7.plot(x_pos, df_sorted['Precision'], 'o-', linewidth=3, markersize=8, 
                label='Precision', color='#FF6B6B', alpha=0.8)
        ax7.plot(x_pos, df_sorted['Recall'], 's-', linewidth=3, markersize=8, 
                label='Recall', color='#4ECDC4', alpha=0.8)
        ax7.plot(x_pos, df_sorted['F1-score'], '^-', linewidth=3, markersize=8, 
                label='F1-Score', color='#45B7D1', alpha=0.8)
        
        ax7.set_title('📊 Metrik Trend Analizi', fontsize=14, fontweight='bold')
        ax7.set_xlabel('Veri Seti (Kalite - Araç)', fontweight='bold')
        ax7.set_ylabel('Metrik Değeri', fontweight='bold')
        ax7.legend(frameon=True, fancybox=True, shadow=True)
        ax7.grid(True, alpha=0.3)
        ax7.set_xticks(x_pos)
        ax7.set_xticklabels([f"{row['Quality']}-{row['Tool']}" for _, row in df_sorted.iterrows()], 
                           rotation=45, ha='right')
        
        return fig
    
    def save_figure(self, fig, filename):
        """Figürü kaydet"""
        output_path = self.output_folder / filename
        fig.savefig(output_path, dpi=300, bbox_inches='tight', 
                   facecolor='white', edgecolor='none', format='png')
        print(f"✅ {filename} başarıyla kaydedildi: {output_path}")
        plt.close(fig)
    
    def generate_all_visualizations(self, csv_filename="metrics.csv"):
        """Tüm görselleştirmeleri oluştur"""
        print("🚀 Görselleştirme işlemi başlatılıyor...")
        
        # Veriyi yükle
        df = self.load_data(csv_filename)
        if df is None:
            return
        
        print(f"📊 Veri özeti:")
        print(f"   • Toplam kayıt: {len(df)}")
        print(f"   • Araçlar: {', '.join(df['Tool'].unique())}")
        print(f"   • Kalite seviyeleri: {', '.join(df['Quality'].unique())}")
        
        # Görselleştirmeleri oluştur
        print("\n🎨 Görselleştirmeler oluşturuluyor...")
        
        # 1. Performans karşılaştırması
        fig1 = self.create_performance_comparison(df)
        self.save_figure(fig1, "01_performance_comparison.png")
        
        # 2. Kalite analizi
        fig2 = self.create_quality_analysis(df)
        self.save_figure(fig2, "02_quality_analysis.png")
        
        # 3. Detaylı metrik tablosu
        fig3 = self.create_detailed_metrics_table(df)
        self.save_figure(fig3, "03_detailed_metrics_table.png")
        
        # 4. Özet dashboard
        fig4 = self.create_summary_dashboard(df)
        self.save_figure(fig4, "04_summary_dashboard.png")
        
        print(f"\n🎉 Tüm görselleştirmeler başarıyla oluşturuldu!")
        print(f"📁 Çıktı klasörü: {self.output_folder.absolute()}")
        
        # Basit rapor
        print(f"\n📋 Hızlı Rapor:")
        print(f"   🏆 En yüksek F1-Score: {df['F1-score'].max():.4f}")
        print(f"   📈 Ortalama Precision: {df['Precision'].mean():.4f}")
        print(f"   🎯 Ortalama Recall: {df['Recall'].mean():.4f}")
        best_overall = df.loc[df['F1-score'].idxmax()]
        print(f"   ⭐ En iyi genel performans: {best_overall['Tool']} ({best_overall['Quality']} kalite)")

def main():
    """Ana fonksiyon"""
    print("🔬 Metrik Görselleştirme Aracı v2.0")
    print("=" * 50)
    
    # Görselleştirici oluştur
    visualizer = MetricsVisualizer(input_folder="input", output_folder="output")
    
    # Tüm görselleştirmeleri oluştur
    visualizer.generate_all_visualizations("metrics.csv")
    
    print("\n✨ İşlem tamamlandı! Çıktı klasörünü kontrol edin.")

if __name__ == "__main__":
    main()