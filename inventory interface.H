// InventoryInterface.h
#pragma once
#include "UObject/Interface.h"
#include "InventoryInterface.generated.h"

UINTERFACE(Blueprintable)
class UInventoryInterface : public UInterface
{
    GENERATED_BODY()
};

class IInventoryInterface
{
    GENERATED_BODY()

public:
    UFUNCTION(BlueprintNativeEvent, BlueprintCallable)
    bool AddItem(TScriptInterface<IItem> Item, int32 Quantity);

    UFUNCTION(BlueprintNativeEvent, BlueprintCallable)
    bool RemoveItem(FName ItemId, int32 Quantity);

    UFUNCTION(BlueprintNativeEvent, BlueprintCallable)
    int32 GetQuantity(FName ItemId) const;
};
