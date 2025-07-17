import pandas as pd
import numpy as np
from statsmodels.api import OLS, add_constant

def get_cross_cor(up_series, down_series, max_lag=10, ini=1):
    """Evalúa la correlación cruzada entre una señal de aguas arriba (up_serie) y una señal de aguas abajo (down_serie)
    
    Args:
        up_serie :  
            serie regularizada 'aguas arriba'
        down_serie :  
            serie regularizada 'aguas abajo' 
        ini : 
            retardo inicial
        max_lag: 
            máximo retardo
    Returns:
        Devuelve un dataframe con retardos (índice) y correlaciones (valor)
    """
    r = []
    for i in range(ini, max_lag + 1):
        u_lagged = up_series.shift(i)
        df = pd.concat([u_lagged, down_series], axis=1).dropna()
        correlation = df.corr().iloc[1, 0]  # Correlation between lagged up_series and down_series
        r.append(correlation)
    return pd.DataFrame(r,index=range(ini, max_lag + 1))

def get_response_time(up_serie : Union[pd.Series,pd.DataFrame], down_serie : Union[pd.Series,pd.DataFrame], max_lag : int =10, ini : int =1) -> int:
    """ Evalúa la correlación cruzada entre una señal de aguas arriba (up_serie) y una señal de aguas abajo (down_serie) y devuelve el índice de retardo que maximiza la correlación para la resolución de las señales consideradas
        Args:
        up_serie :  
            serie regularizada 'aguas arriba'
        down_serie :  
            serie regularizada 'aguas abajo' 
        ini : 
            retardo inicial
        max_lag: 
            máximo retardo
        Returns:
            Devuelve el retraso (índice) para el cual se maximiza la correlación cruzada
    """
    r=get_cross_cor(up_serie,down_serie,max_lag,ini).idxmax()
    return(int(r[0]))

def get_lag_and_linear_fit(up_series : Union[pd.Series,pd.DataFrame], down_series : Union[pd.Series,pd.DataFrame], max_lag : int =10,ini : int =1):
    """Evalúa la correlación cruzada entre una señal de aguas arriba (up_serie) y una señal de aguas abajo (down_serie), determina el tiempo de traslado màs adecuado para la hipótesis de asociación lineal y devuelve el modelo lineal correspondiente
    Args:
        up_serie :  
            serie regularizada 'aguas arriba'
        down_serie :  
            serie regularizada 'aguas abajo' 
        ini : 
            retardo inicial
        max_lag: 
            máximo retardo
        Returns:
            Devuelve el modelo lineal para el cual se maximiza la correlación cruzada
    """
    best_lag = get_response_time(up_series, down_series, max_lag,ini)     
    up_series_shifted = up_series.shift(best_lag)
    df = pd.concat([up_series_shifted, down_series], axis=1).dropna()
    df.columns = ['up_series', 'down_series']
    X = add_constant(df['up_series'])
    y = df['down_series']
    model = OLS(y, X).fit()
    return model

if __name__ == "__main__":
    import sys
