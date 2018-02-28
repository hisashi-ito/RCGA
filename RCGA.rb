#! /usr/bin/env ruby
# -*- coding: utf-8 -*-
#
# 【RCGA】
#
#  概要: Real-coded Genetic Algorithm(RCGA) 計算コード
#
#        パラメータファイルを読み込みRCGAを
#        を実行し試行関数の最適になるパラメータを探索する
#
#  usage: RCGA.rb -c <conf>
#                 -i <parameter file>
#                 -o <output>         
#  引数の説明
#    -c <conf> [必須]
#        設定ファイル
#
#    -i <parmeter file> [必須]
#        パラメータファイル
#
#    -o <output> [必須]
#        結果ファイル  
#
# 遺伝的アルゴリズムの処理の流れ
#
# 1) 初期集団の生成(Initialization)
#    初期集団をランダム変数により生成する
#
# 2) 適応度の評価(Evaluation)
#    試行関数を計算し、条件に満足しているか評価する
#
# 3) 選択(Selection)
#    評価の高い親を選択する
#
# 4) 交配(Crossover)
#    親に似ているデータを生成する
#
# 5) 突然変異(Mutation)
#    進化のために突然変異のデータを少量生成する
#
# 6) 条件判定(Terminate Check)
# 6-0) 適応度の評価(Evaluation)
# 6-1) 収束条件を満足する場合,処理終了
# 6-2) 収束条件を満足しない場合,(2)へ戻る 
#
# 2014.01.23 コードは書き終わり
# 2014.01.26 各種デバッグ作業...
# 2014.01.28 子遺伝子数を制限する
# 2014.01.30 評価値の追加(適合率/再現率/F値)
# 2014.01.30 収束判定条件の変更
#
require 'optparse'
require 'yaml'
require 'logger'

#= 大域定数設定
module CONST
  # 途中経過表示STEP
  STEP=1
end # CONST

#= メッセージ
module MSG
  HELP_MSG =
" usage: RCGA.rb
    -c <conf>  [必須]
        設定ファイル
    -i <input> [必須]
        入力ファイル
　      計算済みの素性値が記載されたtsvファイル
    -o <output> [必須]
        出力ファイル
        探索パラメーター値
"
end # MSG

#= ユーティリティモジュール
module UTIL
  # seed 設定関数
  def set_seed()
    srand(Time.now.to_i + $$)
  end # set_seed
  
  # 範囲付き 浮動小数点型乱数発生関数
  # 引数: x - 下限
  #       y - 上限
  # 返値: 上記範囲内の double の乱数
  def frand(x,y)
    # raise Error if x > y
    rd = ((y.to_f - x.to_f)*rand()) + x
    return rd
  end # frand
end # Utile

#= RCGA クラス
class RealCodedGeneticAlgorithm
  include CONST
  include UTIL
  
  def initialize(input, output, yaml, logger)
    @input = input;
    @output = output;
    @logger = logger
    @yaml = yaml
    
    # 設定ファイル情報を取得
    #= 入力情報
    @parameter_num = @yaml["input"]["parameter_num"].to_i
    @labels = @yaml["input"]["labels"]
    #= パラメータ
    @parameter_ranges = @yaml["parameter_ranges"]
    #= GA 設定
    #== 初期化
    #   母集団数
    @population_num = @yaml["ga"]["initialization"]["population_num"].to_i
    #== 評価
    #    最適化する評価値
    @opt_metric = @yaml["ga"]["evaluation"]["opt_metric"]
    #    収束判定時に平均操作に入れる個体数
    @average_num = @yaml["ga"]["evaluation"]["average_num"].to_i
    #    収束精度
    @tolerance = @yaml["ga"]["evaluation"]["tolerance"].to_f
    #== 選択設定
    #    親の個体数
    @parent_num = @yaml["ga"]["selection"]["parent_num"].to_i
    #    戦略設定
    @strategy = @yaml["ga"]["selection"]["strategy"]
    #    エリート比率
    @elite_ratio = @yaml["ga"]["selection"]["elite_ratio"].to_f
    #== 交配設定
    #    子の個体数
    @crossover_pairs =@yaml["ga"]["crossover"]["crossover_pairs"].to_i
    #    交配確率(ex. 0.8)
    @crossover_prob = @yaml["ga"]["crossover"]["prob"].to_f
    #    交配方法
    @crossover_type = @yaml["ga"]["crossover"]["type"]
    #    交配方法別パラメータ
    #    BLX-alpha の alpha (ex.0.3)
    @alpha = @yaml["ga"]["crossover"]["alpha"].to_f
    
    #== 突然変異
    #    突然変異率(ex. 0.01)
    @mutation_ratio = @yaml["ga"]["mutation"]["prob"].to_f
    #    変異種別設定
    @mutation_type = @yaml["ga"]["mutation"]["type"]
    
    #== 世代
    @max_generation = @yaml["ga"]["max_generation"].to_i
    
    #== 収束判定を実施するために保存する最高精度
    @max_accuracies = []
    
    # 設定チェック
    unless config_chekker()
      # 設定チェックを書いていく
      @logger.fatal("設定に異常があります")
      return 1
    end
    
    #= インタンス変数
    # 正解データ(素性ファイル)
    @relevance = []
    
    # 0) 乱数seed 設定
    set_seed()
  end # initialize
  
  #= 実行関数
  def perform()
    # (0) 計算開始
    #show_infos()

    # (1) 入力データを読み込む
    f = open(@input,"r")
    while line = f.gets
      line.chomp!
      elems = line.split("\t",-1)
      # 0 番目に正解データ
      rel = elems[0]
      # ラベルが不正(未知)の場合は next
      next unless @labels.include?(rel)
      
      # 素性値リスト
      features = elems[1..-1]
      # パラメータ数と素性値数が不一致の場合は next
      next if features.size != @parameter_num
      
      # 正解データに以下の形式でストアする
      # [ 正解ラベル,[ 素性値1,素性値2,...,素性値N ]]
      @relevance.push([rel,features])
    end
    
    # (2) 初期集団作成
    populations = initialization()
    
    # ==> loop
    gen = 0;cflg = false; tol = 0.0; ret = nil;
    while gen < @max_generation
      
      # (3) 評価 
      score_populations = evaluation(populations)
      
      # (3-1) 収束判定
      cflg,tol = convergence?(score_populations)
      if cflg
        # 収束成功: 最終収束パラメータを保存(一番精度の高い結果とパラメータを保存)
        ret =  Marshal.load(Marshal.dump(score_populations[0]))
        break;
      end
      
      # (4) 選択
      selected_individuals = selection(score_populations)
      
      # (5) 交配
      children = crossover(selected_individuals)
      
      # (6) 突然変異
      mutation_children = mutation(children)
      # (6-1) 世代を交代させる
      populations = []
      populations = Marshal.load(Marshal.dump(mutation_children))
      
      # (7) 途中経過出力 (stdout)
      if (gen % STEP) == 0
        if tol.nil?
          @logger.info("first generation")
        else
          @logger.info("generation step: #{gen} tolerance: #{tol}")
        end
      end
      # 世代数(n+1)を進める
      gen += 1
    end
    # <== loop
        
    # (8) 結果表示
    if cflg
      @logger.info("convergence!!")
      @logger.info("max #{@opt_metric}: #{ret[0]}")
      # 収束済み
      # 結果を出力する
      show(ret);
    else
      # 設定世代数を超えても収束しなかった
      @logger.info("設定世代数を超えても収束しませんでした 設定収束精度: #{@tolerance} 最終収束精度: #{tol}")
    end
  end # perform
  
  #= 初期集団生成
  def initialization()
    populations = []
    # 指定された母集団数に従って母集団を作成する
    @population_num.to_i.times{|i|
      params = []
      @parameter_ranges.each do |ary|
        bottom, upper = ary[0].to_f, ary[1].to_f
        # 初期パラメータ(実数)をセットする
        x = frand(bottom, upper)
        params.push(x)
      end
      populations.push(params)
    }
    return populations
  end # initialization
  private :initialization
  
  #= 適応度の評価
  def evaluation(populations)
    # 返却値はスコア付された母集団
    score_populations = []
    # tmp 配列
    tmp = []
    # スコア
    score = 0.0
    # パラメータループ
    populations.each do |params|
      # 指定タイプ毎に計算する評価値を変更する
      case @opt_metric
      # 適合率(平均)
      when "precision" then
        precision_hash = calc_precision(params)
        score = 0.0
        precision_hash.each do |label,ary|
          score += ary[0].to_f / ary[1].to_f
        end
        score = score / precision_hash.size.to_f

      # 再現率(平均)
      when "recall" then
        recall_hash = calc_recall(params)
        score = 0.0
        recall_hash.each do |label,ary|
          score += ary[0].to_f / ary[1].to_f
        end
        score = score / recall_hash.size.to_f

      # F値(平均)
      when "f_mesure" then
        fmesure_hash = calc_fmesure(params)
        score = 0.0
        fmesure_hash.each do |label,f|
          score += f.to_f
        end
        score = score / fmesure_hash.size.to_f
      end

      # [ スコア, [c1,c2,c3,..,cN] ]
      tmp.push([score, params])
    end
    # 評価値毎に昇順にソートする破壊的操作
    tmp.sort!{|a, b| -1*(a[0].to_f <=> b[0].to_f) }
    # 最後に評価値でソートして設定された母集団数にする
    tmp.each_with_index do |ary, idx|
      score_populations.push(ary)
      break if idx >= (@population_num - 1)
    end
    return score_populations
  end # evaluation
  private :evaluation
  
  #== 適合率(precision)
  #   適合率計算用hashを返却する
  def calc_precision(params)
    precision = {}
    # 正解データでループを回す
    @relevance.each do |rel_feat|
      score = 0.0
      rel_label = rel_feat[0]
      # 識別関数計算
      rel_feat[1].each_with_index do |c, idx|
        score += (params[idx].to_f * c.to_f)
      end
      # 識別関数判定
      if score > 0.0
        est_label = @labels[0]
      else
        est_label = @labels[1]
      end
      
      # ラベル毎に以下の情報を保存
      # precision[ 推定ラベル ] = [ 正解数, 推定ラベル総数]
      unless precision.has_key?(est_label)
        # 未登録の場合は初期化
        precision[est_label] = [0, 0]
        # 推定総数をインクリメント
        precision[est_label][1] += 1
        # 正解していたら正解数をインクリメント
        if rel_label == est_label
          precision[est_label][0] += 1
        end
      else
        # 推定総数をインクリメント
        precision[est_label][1] += 1
        # 正解していたら正解数をインクリメント
        if rel_label == est_label
          precision[est_label][0] += 1
        end
      end
    end # relevance
    return precision
  end # calc_precision
  private :calc_precision
    
  #== 再現率(recall)
  #   再現率計算用hashを返却する
  def calc_recall(params)
    recall = {}
    # 正解データでループを回す
    @relevance.each do |rel_feat|
      score = 0.0
      rel_label = rel_feat[0]
      # 評価関数計算
      rel_feat[1].each_with_index do |c, idx|
        score += (params[idx].to_f * c.to_f)
      end
      # 識別関数判定
      if score > 0.0
        est_label = @labels[0]
      else
        est_label = @labels[1]
      end

      # ラベル毎に以下の情報を保存
      # recall[ 正解ラベル ] = [ 正解数, 正解ラベル総数]
      unless recall.has_key?(rel_label)
        # 未登録の場合は初期化
        recall[rel_label] = [0, 0]
        # 推定総数をインクリメント
        recall[rel_label][1] += 1
        # 正解していたら正解数をインクリメント
        if rel_label == est_label
          recall[rel_label][0] += 1
        end
      else
        # 推定総数をインクリメント
        recall[rel_label][1] += 1
        # 正解していたら正解数をインクリメント
        if rel_label == est_label
          recall[rel_label][0] += 1
        end
      end
    end # relevance
    return recall
  end # calc_recall
  private :calc_recall
  
  #== F値(F-mesure)
  #   F値計算用hashを返却する
  def calc_fmesure(params)
    fmesure = {}
    # 適合率計算(precision)
    precision = calc_precision(params)
    
    # 再現率計算(recall)
    recall = calc_recall(params)
    
    # ラベル毎にF値 を計算する
    @labels.each do |label|
      p = precision[label][0].to_f / precision[label][1].to_f
      r = recall[label][0].to_f / recall[label][1].to_f
      f =  (2.0 * p * r) / ( p + r )
      fmesure[label] = f
    end
    return fmesure
  end # calc_fmesure
  private :calc_fmesure
    
  #== 収束判定
  #   収束判定条件を以下のように定義する
  #
  #   母集団の内で一定数のスコア上位個体の中で評価値の平均を計算し
  #   平均評価値の値の変化比率がある一定閾値以下になった場合に、
  #   収束判定したと見なす
  #
  def convergence?(score_populations)
    # 最高スコアを取得
    max_score = score_populations[0][0].to_f
    # 平均スコア
    average_score = 0.0;
    score_populations.each_with_index do |sp,idx|
      average_score += sp[0].to_f
      # スコア上位の指定された母集団数でしか評価しない
      break if idx >= @average_num
    end
    
    # 初回時の処理
    if @max_accuracies[-1].nil?
      # 平均スコアをインスタンス変数へ保存
      @max_accuracies.push(average_score)
      return [false, nil]
    end
    # 一つ前の平均スコア値を取得する
    pre_average_score = @max_accuracies[-1]
    
    # 差分を計算し、収束判定を実施する
    tol = ((average_score - pre_average_score).abs / pre_average_score)
    
    # カレントスコアを最終行に追加
    @max_accuracies.push(average_score)

    if tol < @tolerance
      # 収束成功!
      # Boolean値とtorelanceを返却
      return [true, tol]
    else
      # 収束まだ...
      # 収束しなくてもtorelanceを返却
      return [false, tol]
    end
  end # convergence?
  private :convergence?
  
  #= 選択
  def selection(score_populations)
    # 選択された個体集合
    selected_individuals = []
    
    # 戦略毎に母集団の中から親を選択する
    case @strategy
      
    # [エリート戦略]
    when 'elite' then
      # 評価値(スコア)の高い順番に規定数抽出する
      score_populations.each_with_index do |ary, idx|
        selected_individuals.push(ary)
        break if idx >= (@parent_num - 1)
      end

    # [ルーレット戦略]
    when 'roulette' then
      # (0) Temp に母集団をコピーする
      tmp = Marshal.load(Marshal.dump(score_populations))
      # (1) 評価値の総和を求める
      sum = 0.0;
      score_populations.each do |ary|
        sum += ary[0].to_f
      end
      # (2) (0..sum) の間の数をランダムに選択する
      rs = frand(0.0, sum)
      # (3) 個体をランダムに選択する
      th_score = 0.0;
      loop {
        r_idx = rand(tmp.size)
        th_score += tmp[r_idx][0]
        # 可算スコアが指定スコア以上であれば,break
        break if th_score >= rs
        # 可算スコアが指定スコア以下であれば情報を保存
        selected_individuals.push(tmp[r_idx])
        tmp.delete_at(r_idx)
        # 親世代が指定数以上に達した場合はbreak
        break if selected_individuals.size > @parent_num
      }

    # [ハイブリッド戦略]
    #  エリート戦略で規定数を取得後、
    #  残りの個体に対してルーレット選抜を実施
    when 'hybrid' then
      # (0) Temp に母集団をコピーする
      tmp = Marshal.load(Marshal.dump(score_populations))
      # (1) 最初にエリート戦略でどれぐらい選択するか計算する
      #     設定で指定された率を親抽出数に掛けたものとする
      elite_num = (@parent_num * @elite_ratio).to_i
      tmp.each_with_index do |ary, idx|
        selected_individuals.push(ary)
        break if idx >= (elite_num - 1)
      end
      # (2) エリート戦略数で該当しなかったもんからルーレット戦略
      #     で選択する
      tmp2 = Marshal.load(Marshal.dump(tmp[elite_num..-1]))
      # (3) 評価値の総和を求める
      sum = 0.0;
      tmp2.each do |ary|
        sum += ary[0].to_f
      end
      # (4) (0..sum) の間の数をランダムに選択する
      rs = frand(0.0, sum)
      # (5) 個体をランダムに選択する
      th_score = 0.0;
      loop {
        r_idx = rand(tmp2.size)
        th_score += tmp2[r_idx][0]
        # スコアが指定スコア以上であれば,break
        break if th_score >= rs
        # スコアが指定スコア以下であれば情報を保存
        selected_individuals.push(tmp2[r_idx])
        tmp2.delete_at(r_idx)
        # 親世代が指定数以上に達した場合はbreak
        break if selected_individuals.size > @parent_num
      }
    else
      # (ここには入らない）
    end
    return selected_individuals
  end # selection
  private :selection
  
  #= 交配
  def crossover(selected_individuals)
    # 次世代個体(子供たち)
    children = []
    
    # 交配方法毎に処理を変更する
    case @crossover_type

    # [BLX-alpha交配]
    when "BLX-alpha" then
      children = blx_alpha(selected_individuals)
      
    # [シンプレックス交配]
    when "SPX" then
      children = simplex(selected_individuals)
    else
      # no effect
    end
    return children
  end # crossover
  private :crossover
  
  #== [BLX-alpha 交配]
  def blx_alpha(selected_individuals)
    parents = [];
    children = [];
    indi_size = selected_individuals.size
    
    # selected_individuals の中からランダムに 
    # crossover_pairs 個の親のセット[a,b]を作る
    loop{
      # 全部のループを回すのではなく交配確率に依存
      next if frand(0.0,1.0) > @crossover_prob
      # ランダムで親(a),(b) の２つを選択
      rn_a = rand(indi_size)
      rn_b = nil
      # ホモ交配は許容しない(双方が異なる配列要素ができるまでloop を回す)
      loop{
        rn_b = rand(indi_size)
        break if rn_b != rn_a
      }
      # パラメータを保存
      # [accurary,[c1,c2,c3...cN]]
      parent_a = selected_individuals[rn_a[1]]
      parent_b = selected_individuals[rn_b[1]]
      parents.push([parent_a, parent_b])
      break if parents.size >= @crossover_pairs
    }
    
    # 交配する親の個体について2 組のすべてのペアを作成する
    parents.each{|ary|
      # A x B のペアを作成する
      # A = [ca1,ca2,...caN]
      # B = [cb1,cb2,...cbN]
      parent_a = ary[0][1]
      parent_b = ary[1][1]
      
      # 新規に生成する子のパラメータの初期化
      child = []
      
      # 2の親のパラメータ毎にループ
      parent_a.each_with_index do |pa, idx|
        pb = parent_b[idx]
        # BLX-alphaの範囲(下限,上限)を計算
        bottom, upper = blx_alpha_range(pa,pb)
        # 子のパラメータを生成
        c = frand(bottom, upper)
        child.push(c)
      end
      # 子のパラメータがすべて計算されたので保存する
      children.push(child)
    }
    
    # 最後に元の母集合を子集合の配列に結合する
    # 結合後の母集団のサイズは 
    # population_num + crossover_pairs(1:1生成なら)
    selected_individuals.each do |ary|
      # ary[1]にパラメータが保持されている
      children.push(ary[1])
    end
    return children
  end
  private :blx_alpha
  
  #=== [BLX-alpha 交配 範囲計算]
  def blx_alpha_range(x1,x2)
    dx = (x1 - x2).abs
    min_x = [x1, x2].min
    max_x = [x1, x2].max
    # 次世代範囲x軸
    min_cx = min_x - (@alpha * dx)
    max_cx = max_x + (@alpha * dx)
    return [min_cx, max_cx]
  end # blx_alpha_range
  
  #== [シンプレックス交配]
  def simplex(selected_individuals)
    # epsilon = sqrt(n+2)
    epslion = Math.sqrt(@parameter_num + 2.0)
    parents = [];
    children = [];
    
    # 親を規定する(交配確率を計算し一定値以上を親とする)
    selected_individuals.each do |ary|
      # 交配確率で指定値以上の値であれば交配を実施しない
      next if frand(0.0,1.0) > @crossover_prob
      # [accuracy, [c1,c2,c3,...,cN]]
      parents.push(ary[1])
    end
    
    # シンプレックス交配は遺伝子長(n)の場合
    # n+1 個の親個体により1個の個を生成する
    # ary = [p1,p2,p3,,,pN+1]
    parents.combination(@parameter_num.to_i + 1){|ary|
      g_vec = [];s_vec = []
      # 重心ベクトル計算
      ary.each_with_index do |params,idx|
        @parameter_num.times do |jdx|
          g_vec[jdx] = 0.0 if g_vec[jdx].nil?
          g_vec[jdx] += params[jdx]
        end
      end
      g_vec.map!{|elem| elem/(@parameter_num + 1).to_f}
      
      # S ベクトル計算
      ary.each_with_index do |params, idx|
        s_vec[idx] = []
        params.each_with_index do |c, jdx|
          s_vec[idx][jdx] = g_vec[jdx] + (epslion * (c - g_vec[jdx]))
        end
      end
      
      # 子パラメータ生成
      cpara = []
      @parameter_num.times{|idx|
        if idx == 0
          # Ci = Si  (i=0)
          cpara.push(s_vec[idx])
        else
          # Ci+1 = Si+1 + ri*(Si+1-Si+Ci))
          # r = ( u[0,1] )^1/i
          r = frand(0.0,1.0)**(1.0/idx.to_f)
          
          # デバッグ中
          
          c = s_vec[idx] + (r * (s_vec[idx] - s_vec[idx-1] + cpara[-1]))
          cpara << c
        end
        # 最後に子パラメータを配列
        children.push(cpara)
      }
    }
    return children    
  end # simplex
  private :simplex
  
  #= 突然変異
  #  次世代パラメータに対してある一定の確率で突然変異を発生させる
  #  本関数の中で２つの突然変異型の挙動を実装する
  #  ・一様突然変異
  #  ・境界突然変異
  #
  def mutation(children)
    mutation_children = []
    children.each do |param|
      # param = [c1,c2,c3...cN]
      if frand(0.0, 1.0) > @crossover_prob
        # [0..1] 変数で突然変異率よりも大きい場合は突然変異しない
        mutation_children.push(param)
      else
        # 突然変異を実施
        # 何箇所突然変異させるかは乱数で決める
        m_param = Marshal.load(Marshal.dump(param))
        rand(@parameter_num + 1).times do |idx|
          # どこの遺伝子を突然変異させるかも乱数で決める
          mdx = rand(@parameter_num)
          bottom, upper = @parameter_ranges[mdx]
          
          # 突然変異タイプ毎の処理
          # [一様突然変異]
          #  範囲内を乱数値で新パラメータを設定する
          if @mutation_type == 'uniform'
            m_param[mdx] = frand(bottom, upper)

          # [境界突然変異]
          # 境界値を設定する、上限下限のどちらを設定するかは乱数で決める
          elsif @mutation_type == 'boundary'
            if frand(0.0, 1.0) > 0.5 # 確率はeven
              m_param[mdx] = upper.to_f
            else
              m_param[mdx] = bottom.to_f
            end
            # ここは入らない
          end
        end
        # 突然変異済みにパラメータに追加
        mutation_children.push(m_param)
      end
    end
    return mutation_children
  end # mutation
  private :mutation
  
  #= show_infos
  #  計算前に各種情報を出す
  def show_infos()
    @logger.info("/// CONF SETTEINGS ///")
    @logger.info("パラメータ数: #{@parameter_num}")
    @logger.info("指定ラベル: #{@labels.join("\s")}")
    @logger.info("探索用パラメータ設定:")
    @parameter_ranges.each do |ary|
      bottom, upper = ary[0],ary[1]
      @logger.info("#{bottom},#{upper}")
    end
    @logger.info("//// GA /////")
    @logger.info("最大世代数: #{@max_generation}")
    @logger.info("母集団数数: #{@population_num}")
    @logger.info("収束精度: #{@tolerance}")
  end # showw_infos
  private :show_infos
  
  #= show
  #  最終出力
  def show(ret)
    accuracy = ret[0]
    params = ret[1]
    o = open(@output,"w")
    o.puts "最高精度: #{accuracy}"
    o.puts "*** parameters ***"
    params.each_with_index do |c, idx|
      o.puts "c_#{idx}: #{c}"
    end
    o.close
  end # show
  private :show
  
  #= config cheker
  def config_chekker()
    # 母集団が親の個体数よりも小さい場合
    if @population_num < @parent_num
      @logger.fatal("母集団数が親の個体数よりも小さいです 母体数:#{@population_num} 親個体数:#{@parent_num}")
      return false
    end
    # 指定パラメーター数とパラメータのレンジ設定が異なります
    if @parameter_ranges.size != @parameter_num
      @logger.fatal("指定パラメータ数とパラメータ範囲数が異なります パラメータ範囲数: #{parameter_ranges.size} パラメータ数: #{parameter_num}")
      return false
    end
    
    # 全部の設置がO.K.であればtrurを返却
    return true
  end # config_chekker
  private :config_chekker
end # real coded genetic algorithm


# メイン関数
def main(argv)
  include MSG
  input = nil;output = nil;conf = nil;hflg = false
  # option parser
  OptionParser.new {|opt|
    opt.on('-i input'){|v| input = v}
    opt.on('-o output'){|v| output = v}
    opt.on('-c conf'){|v| conf = v}
    opt.on('-h'){|v| hflg = v}
    opt.parse!(argv)
  }
  
  # ヘルプ表示
  if hflg
    puts HELP_MSG;
    return 1
  end
  
  # 引数チェック
  if input.nil? || output.nil? || conf.nil?
    $stderr.puts HELP_MSG;
    return 1
  end
  
  # 設定ファイルを読み込む
  yaml = YAML.load_file conf
  # logger をセット
  logger = Logger.new(STDERR)
  logger_level = yaml["log"]["level"].to_s
  if logger_level == 'INFO'
    logger.level = Logger::INFO
  elsif logger_level == 'WARN'
    logger.level = Logger::WARN
  elsif logger_level == 'ERROR'
    logger.level = Logger::ERROR
  elsif logger_level == 'FATAL'
    logger.level = Logger::FATAL
  else
    $stderr.puts "error log level"
    return 1
  end
  # 処理開始
  logger.info("*** start RCGA ***")
  
  # Real-coded GA インスタンス作成
  ga = RealCodedGeneticAlgorithm.new(input, output, yaml, logger);
  
  # 探索開始
  ga.perform();
  
  logger.info("*** stop RCGA ***")
  return 0
end

if __FILE__ == $0
  main(ARGV)
end
