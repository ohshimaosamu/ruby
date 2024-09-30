def compute_Lb_t(la, lb, a, b, d)
  lb_t = Array.new(b){Array.new(b,0)} # 輝度配列Lb_tを初期化
  (0...b).each do |y|  
    (0...b).each do |x|
      xa=1-(d-x)
      if xa >= (b-a)/2.0 && xa < a+(b-a)/2.0 && y>=(b-a)/2.0 && y < a+(b-a)/2.0
        lb_t[x][y] =  lb[x][y] * la[xa-(b-a)/2.0][y-(b-a)/2.0]
      else
        lb_t[x][y] = lb[x][y]
      end
    end
  end
  return lb_t
end
def rotate_point(x, y, theta)
  # 座標を極座標に変換
  r = Math.sqrt(x**2 + y**2)
  phi = Math.atan2(y, x)
  # 回転後の角度を計算
  phi_prime = phi + theta
  # 極座標を直交座標に変換
  x_prime = r * Math.cos(phi_prime)
  y_prime = r * Math.sin(phi_prime)
  return x_prime, y_prime
end
def point_inside_ellipse?(x, y, a, b) #点(x, y)が楕円の内部にあるかどうかを判定  
  return ((x**2) / (a**2) + (y**2) / (b**2)) < 1
end

def get_la(a,a1,b1,theta)
	rad = -1*Math::PI * theta/180.0  # 度回転
	la=Array.new(a){Array.new(a,0)}
	(0...a).each do |x|
 	 	(0...a).reverse_each do |y|
    		xd, yd = rotate_point((x-a/2.0).to_f, (y-a/2.0).to_f, rad)
    		if point_inside_ellipse?(xd, yd, a1/2.0, b1/2.0)
    	  		la[x][y]=0  #楕円の内部にあります
    		else
      			la[x][y]=1  #楕円の外部にあります
    		end
  		end
	end
  return la
end

def ld(mu,bnd)  #bnd=0,1,2,3 for UBVRI 近似式によるLDの計算　半径に対応するmu　バンド名
#Claret 2022 logg=0.0 Teff=3750 Z=-0.2 Vt=2.0 ATLAS
g=[0.8379, 0.7975, 0.8112, 0.8235]
h=[1.0344, 1.4398, 1.2091, 0.9642]
	return 1-g[bnd]*(1 - mu**h[bnd])
end

def lumi(a,band) # bandによる移動量aでの食の光量
  la = Array.new(a){Array.new(a,0)}
  (0...a).each do |x|
    (0...a).each do |y|
      r=Math.sqrt((x-a/2.0)**2+(y-a/2.0)**2)
      if r<(a-2)/2.0 
        la[x][y]=ld(1-Math.sqrt((x-a/2.0)**2+(y-a/2.0)**2.0)/a/2.0, band)
#        la[x][y]=1 
      else
        la[x][y]=0
      end
    end
  end
  return la
end

#main routine
# 例: 正方形a, bの輝度配列とそれぞれの1辺のピクセル数を与える
a = 240		#a　小惑星マップ範囲の1辺の長さ
a1 = 200	#a1 長軸長not半長軸
b1 = 185	#b1 短軸長not半短軸
theta = 30	#theta = 60  #楕円の長軸がx軸から時計回りにtheta度回転しているとする
b = 220  #radus of Betelgeuse
# 正方形aと正方形bの輝度配列la(x,y)とlb(x,y)
 la = get_la(a,a1,b1,theta)
lb=[]
lb0=[]
(0..3).each do |band|
 lb[band]=lumi(b,band)
# 正方形bの輝度配列Lb_tを計算
  lb0[band]=compute_Lb_t(la, lb[band], a, b, (a+b)/2.0).flatten.sum
end
  min=1.0
print "#a=",a," a1=",a1," b1=",b1," theta=",theta," b=",b,"\n"
puts "#a:size of asteroide mask, a1:dia of major axis, b1:dia of minor axis, theta:angle of major axis from fowerding direction, b:size of stellar map"
puts "#d    	U           	B            	V            	R"
(-(a+b)/2..(a+b)/2+1).each do |d|
  printf "%2.3f ", 2*d.to_f/(a+b).to_f
  (0..3).each do |band|
    lb_t = compute_Lb_t(la, lb[band], a, b, d)
    l_total=lb_t.flatten.sum/lb0[band]
    print l_total," "
  end
  print "\n"
end

