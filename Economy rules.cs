public class ExchangeRate
{
    public string FromItemGlobalId { get; set; }
    public string ToItemGlobalId { get; set; }
    public double Rate { get; set; } // 1 From = Rate To
}

public class CrossGameEconomy
{
    public Dictionary<string, ExchangeRate> Rates { get; set; } = new();

    public bool TryConvert(string fromGlobalId, string toGlobalId, long fromQuantity, out long toQuantity)
    {
        toQuantity = 0;
        string key = $"{fromGlobalId}->{toGlobalId}";
        if (!Rates.TryGetValue(key, out var rate))
            return false;

        toQuantity = (long)(fromQuantity * rate.Rate);
        return true;
    }
}
