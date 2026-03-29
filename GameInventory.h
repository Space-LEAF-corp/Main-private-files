// GameInventory.h
#pragma once
#include "CoreMinimal.h"
#include "InventoryInterface.h"
#include "GameInventory.generated.h"

UCLASS(Blueprintable)
class UGameInventory : public UObject, public IInventoryInterface
{
    GENERATED_BODY()

protected:
    UPROPERTY()
    TMap<FName, int32> Items;

public:
    virtual bool AddItem_Implementation(TScriptInterface<IItem> Item, int32 Quantity) override;
    virtual bool RemoveItem_Implementation(FName ItemId, int32 Quantity) override;
    virtual int32 GetQuantity_Implementation(FName ItemId) const override;
};
